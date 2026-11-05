#!/usr/bin/env python3
"""Build a censored two-factor probabilistic binding representation.

The observation model is y_ij ~ N(mu_j + W_j z_i, sigma_j^2), z_i ~ N(0, I_2),
for j in {Ki, Kd, IC50}. Inequality-qualified concentrations are represented as
one-sided bounds after the -log10(molar) transformation. Posterior means of z_i
form the two-dimensional drug-target representation.
"""

from __future__ import annotations

import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import log_ndtr, ndtr


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from drug_synergy.paths import ANNOTATIONS, INPUTS, RESULTS  # noqa: E402

RAW_INPUT = INPUTS / "binding_assays_raw_harmonized.tsv"
ALMANAC = ANNOTATIONS / "almanac_drugs.csv"
OUTPUT = RESULTS / "pfa_binding"
ASSAYS = ["Ki", "Kd", "IC50"]
RANK = 2


def target_name(value):
    return " ".join(str(value).strip().split())


def concentration_to_p(value, units="nM"):
    try:
        number = float(str(value).replace(",", "").strip())
    except (TypeError, ValueError):
        return None
    factor = {
        "pm": 1e-12, "nm": 1e-9, "um": 1e-6, "µm": 1e-6,
        "mm": 1e-3, "m": 1.0,
    }.get(str(units).strip().lower())
    if factor is None or not math.isfinite(number) or number <= 0:
        return None
    return -math.log10(number * factor)


def parse_text_measurement(value):
    text = "" if pd.isna(value) else str(value).strip()
    match = re.match(r"^\s*([<>]=?|=|~)?\s*([0-9,.]+(?:e[-+]?\d+)?)", text, re.I)
    if not match:
        return None, None
    relation = match.group(1) or "="
    return relation, concentration_to_p(match.group(2), "nM")


def p_kind(relation):
    relation = str(relation or "=").strip()
    if relation.startswith(">"):
        return "upper"  # concentration > c implies p < p(c)
    if relation.startswith("<"):
        return "lower"  # concentration < c implies p > p(c)
    return "exact"


def add(rows, drug, target, assay, value, relation, units, source):
    if assay not in ASSAYS or pd.isna(target) or str(target).strip().lower() == "nan":
        return
    p = concentration_to_p(value, units)
    if p is None:
        return
    rows.append({
        "almanac_drug": str(drug), "target_name": target_name(target),
        "assay_type": assay, "p_value": p, "relation_original": str(relation or "="),
        "p_bound_type": p_kind(relation), "source": source,
    })


def load_records():
    """Load the frozen harmonized measurements used in the manuscript."""
    raw = pd.read_csv(RAW_INPUT, sep="\t").drop_duplicates([
        "almanac_drug", "target_name", "assay_type", "p_value",
        "relation_original", "source",
    ])
    return raw


def collapse_bounds(raw):
    records = []
    keys = ["almanac_drug", "target_name", "assay_type"]
    for key, group in raw.groupby(keys, sort=False):
        exact = group.loc[group.p_bound_type == "exact", "p_value"].to_numpy(float)
        lower = group.loc[group.p_bound_type == "lower", "p_value"].to_numpy(float)
        upper = group.loc[group.p_bound_type == "upper", "p_value"].to_numpy(float)
        row = dict(zip(keys, key))
        row.update({
            "exact_p": float(np.median(exact)) if len(exact) else np.nan,
            "lower_p": float(np.max(lower)) if len(lower) else -np.inf,
            "upper_p": float(np.min(upper)) if len(upper) else np.inf,
            "n_exact": len(exact), "n_lower_censored": len(lower),
            "n_upper_censored": len(upper), "n_records": len(group),
            "sources": ";".join(sorted(group.source.unique())),
        })
        if len(exact):
            row["observation_type"] = "exact"
        elif len(lower) and len(upper) and row["lower_p"] < row["upper_p"]:
            row["observation_type"] = "interval"
        elif len(lower) and not len(upper):
            row["observation_type"] = "lower"
        elif len(upper) and not len(lower):
            row["observation_type"] = "upper"
        else:
            row["observation_type"] = "inconsistent_bounds"
        records.append(row)
    return pd.DataFrame(records)


def truncated_expectation(mean, sigma, lower, upper):
    """Mean of N(mean,sigma^2) truncated to [lower, upper]."""
    if np.isneginf(lower) and np.isposinf(upper):
        return mean
    if np.isneginf(lower):
        b = (upper - mean) / sigma
        ratio = math.exp(-0.5 * b * b - 0.5 * math.log(2 * math.pi) - log_ndtr(b))
        return mean - sigma * ratio
    if np.isposinf(upper):
        a = (lower - mean) / sigma
        ratio = math.exp(-0.5 * a * a - 0.5 * math.log(2 * math.pi) - log_ndtr(-a))
        return mean + sigma * ratio
    a, b = (lower - mean) / sigma, (upper - mean) / sigma
    za, zb = math.exp(log_ndtr(a)), math.exp(log_ndtr(b))
    denom = max(zb - za, 1e-14)
    phi_a = math.exp(-0.5 * a * a) / math.sqrt(2 * math.pi)
    phi_b = math.exp(-0.5 * b * b) / math.sqrt(2 * math.pi)
    value = mean + sigma * (phi_a - phi_b) / denom
    return float(np.clip(value, lower, upper))


def fit_pfa(exact, lower, upper, mask, assay_sigma, rank=2, max_iter=80, ridge=1e-3):
    n, assays = exact.shape
    pseudo = exact.copy()
    for i in range(n):
        for j in range(assays):
            if not mask[i, j] or np.isfinite(pseudo[i, j]):
                continue
            lo, hi = lower[i, j], upper[i, j]
            if np.isfinite(lo) and np.isfinite(hi):
                pseudo[i, j] = (lo + hi) / 2
            elif np.isfinite(lo):
                pseudo[i, j] = lo + 0.5
            else:
                pseudo[i, j] = hi - 0.5
    mu = np.array([np.nanmedian(pseudo[:, j]) for j in range(assays)])
    centered = np.where(mask, pseudo - mu, 0.0)
    u, s, _ = np.linalg.svd(centered, full_matrices=False)
    z = u[:, :rank] * s[:rank]
    z /= np.std(z, axis=0, keepdims=True) + 1e-8
    w = np.zeros((assays, rank))
    sigma = np.asarray(assay_sigma, dtype=float).copy()
    previous = np.inf
    posterior_sd = np.zeros((n, rank))
    posterior_cov = np.repeat(np.eye(rank)[None, :, :], n, axis=0)
    for iteration in range(max_iter):
        for j in range(assays):
            indices = np.flatnonzero(mask[:, j])
            zj = z[indices]
            sum_z = zj.sum(axis=0)
            sum_zz = zj.T @ zj + posterior_cov[indices].sum(axis=0)
            normal = np.block([
                [np.array([[len(indices)]], dtype=float), sum_z[None, :]],
                [sum_z[:, None], sum_zz],
            ])
            normal += np.diag([0.0] + [ridge] * rank)
            response = np.concatenate([
                [pseudo[indices, j].sum()],
                zj.T @ pseudo[indices, j],
            ])
            coef = np.linalg.solve(normal, response)
            mu[j], w[j] = coef[0], coef[1:]

        patterns = np.unique(mask, axis=0)
        for pattern in patterns:
            indices = np.flatnonzero(np.all(mask == pattern, axis=1))
            obs = np.flatnonzero(pattern)
            precision = np.eye(rank)
            for j in obs:
                precision += np.outer(w[j], w[j]) / sigma[j] ** 2
            covariance = np.linalg.inv(precision)
            weighted_loadings = w[obs] / sigma[obs, None] ** 2
            rhs = (pseudo[np.ix_(indices, obs)] - mu[obs]) @ weighted_loadings
            z[indices] = rhs @ covariance.T
            posterior_cov[indices] = covariance
            posterior_sd[indices] = np.sqrt(np.diag(covariance))

        prediction = mu + z @ w.T
        for j in range(assays):
            censored = mask[:, j] & ~np.isfinite(exact[:, j])
            for i in np.flatnonzero(censored):
                pseudo[i, j] = truncated_expectation(
                    prediction[i, j], sigma[j], lower[i, j], upper[i, j]
                )
        standardized = (pseudo - prediction) / sigma
        loss = float(np.mean(standardized[mask] ** 2))
        if abs(previous - loss) < 1e-7:
            break
        previous = loss

    # Orthogonal rotations leave the combined downstream cosine similarity unchanged.
    order = np.argsort(-np.linalg.norm(w, axis=0))
    z, posterior_sd, w = z[:, order], posterior_sd[:, order], w[:, order]
    return z, posterior_sd, mu, w, sigma, pseudo, iteration + 1


def scaled_cosine(frame):
    values = frame.fillna(0.0).to_numpy(float)
    norms = np.linalg.norm(values, axis=1)
    denom = np.outer(norms, norms)
    cosine = np.divide(values @ values.T, denom, out=np.zeros_like(denom), where=denom > 0)
    np.fill_diagonal(cosine, 1.0)
    return pd.DataFrame((cosine + 1.0) / 2.0, index=frame.index, columns=frame.index)


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    raw = load_records()
    collapsed = collapse_bounds(raw)
    usable = collapsed[collapsed.observation_type != "inconsistent_bounds"].copy()
    collapsed.to_csv(OUTPUT / "binding_assays_collapsed.tsv", sep="\t", index=False)

    index = usable[["almanac_drug", "target_name"]].drop_duplicates()
    index["row"] = np.arange(len(index))
    row_map = {(r.almanac_drug, r.target_name): int(r.row) for r in index.itertuples()}
    assay_map = {name: j for j, name in enumerate(ASSAYS)}
    shape = (len(index), len(ASSAYS))
    exact = np.full(shape, np.nan)
    lower = np.full(shape, -np.inf)
    upper = np.full(shape, np.inf)
    mask = np.zeros(shape, dtype=bool)
    for row in usable.itertuples(index=False):
        i, j = row_map[(row.almanac_drug, row.target_name)], assay_map[row.assay_type]
        exact[i, j], lower[i, j], upper[i, j] = row.exact_p, row.lower_p, row.upper_p
        mask[i, j] = True

    exact_raw = raw[raw.p_bound_type == "exact"].copy()
    exact_raw["group_median"] = exact_raw.groupby(
        ["almanac_drug", "target_name", "assay_type"]
    ).p_value.transform("median")
    exact_raw["replicate_residual"] = exact_raw.p_value - exact_raw.group_median
    assay_sigma = np.asarray([
        max(float(np.sqrt(np.mean(
            exact_raw.loc[exact_raw.assay_type == assay, "replicate_residual"] ** 2
        ))), 0.25)
        for assay in ASSAYS
    ])
    z, zsd, mu, w, sigma, pseudo, iterations = fit_pfa(
        exact, lower, upper, mask, assay_sigma
    )
    latent = index.drop(columns="row").copy()
    latent["latent_z1"], latent["latent_z2"] = z[:, 0], z[:, 1]
    latent["posterior_sd_z1"], latent["posterior_sd_z2"] = zsd[:, 0], zsd[:, 1]
    latent["observed_assay_count"] = mask.sum(axis=1)
    latent.to_csv(OUTPUT / "pfa2_latent.tsv", sep="\t", index=False)
    parameters = pd.DataFrame({
        "assay_type": ASSAYS, "intercept": mu, "loading_z1": w[:, 0],
        "loading_z2": w[:, 1], "residual_sigma": sigma,
        "effective_precision_weight": 1 / sigma ** 2, "observed_cells": mask.sum(axis=0),
    })
    parameters.to_csv(OUTPUT / "pfa2_parameters.tsv", sep="\t", index=False)

    drugs = pd.read_csv(ALMANAC).drug.astype(str).tolist()
    z1 = latent.pivot(index="almanac_drug", columns="target_name", values="latent_z1").reindex(drugs)
    z2 = latent.pivot(index="almanac_drug", columns="target_name", values="latent_z2").reindex(drugs)
    z1.columns = [f"{column}::z1" for column in z1.columns]
    z2.columns = [f"{column}::z2" for column in z2.columns]
    profile = pd.concat([z1, z2], axis=1)
    profile.to_csv(OUTPUT / "pfa2_drug_target_profile.tsv", sep="\t")
    similarity = scaled_cosine(profile)
    similarity.to_csv(OUTPUT / "pfa2_binding_similarity.csv")
    has_profile = profile.notna().any(axis=1)
    pair_mask = pd.DataFrame(
        np.outer(has_profile.astype(int), has_profile.astype(int)),
        index=profile.index, columns=profile.index,
    )
    np.fill_diagonal(pair_mask.values, 0)
    pair_mask.to_csv(OUTPUT / "pfa2_binding_pair_mask.csv")

    print(f"raw unique measurements: {len(raw):,}")
    print(raw.p_bound_type.value_counts().to_string())
    print(f"collapsed assay cells: {len(collapsed):,}")
    print(collapsed.observation_type.value_counts().to_string())
    print(f"drug-target rows: {len(latent):,}; fitted iterations: {iterations}")
    print(parameters.to_string(index=False))


if __name__ == "__main__":
    main()
