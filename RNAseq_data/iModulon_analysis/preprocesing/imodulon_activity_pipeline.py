#!/usr/bin/env python3
"""Infer PRECISE411 iModulon activities for Figure 6C and Table S9.

Run from the manuscript project root:
    python RNAseq_data/iModulon_analysis/preprocesing/imodulon_activity_pipeline.py \
        --config RNAseq_data/iModulon_analysis/preprocesing/config.json

See README_iModulon.md beside this script for input locations and sources.
Paths in config.json are relative to that configuration file, so no personal
computer paths or change of working directory are needed.

Analysis flow:
1. Read replicate-level log2(TPM + 1), sample groups and the PRECISE411 reference.
2. Match gene IDs and subtract each gene's mean in untreated controls.
3. Infer all 116 activities together: A = pinv(M) @ centered expression.
4. Compare sample activities using group OLS contrasts and two-sided t tests.
5. Apply BH across all 580 tests; call changes at q < 0.05, |delta| > 5,
   with adequate reference coverage. Export the complete result table.

The projection follows pymodulon.util.infer_activities, implemented here with
NumPy; installing pymodulon is not required. Primary inference uses biological
replicates, not PyModulon's historical log-normal reference-null DiMA test.
Original calculations, options and diagnostic outputs are retained. Optional
reference-null and target-exclusion analyses are disabled in the supplied config.
"""
from __future__ import annotations

import argparse
from itertools import combinations
import hashlib
import json
import logging
from pathlib import Path
import re
import sys
import warnings

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import linkage, leaves_list
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

VERSION = "1.0.0"
LOG = logging.getLogger("imodulon")


def read_table(path):
    """CSV or TSV, including gzip; identifiers remain strings."""
    sep = "\t" if ".tsv" in str(path).lower() else ","
    return pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False)


def matrix(path, index_column=None):
    t = read_table(path)
    key = index_column or t.columns[0]
    if key not in t:
        raise ValueError(f"Index column {key!r} absent in {path}")
    t[key] = t[key].str.strip()
    if (t[key] == "").any() or t[key].duplicated().any():
        raise ValueError(f"Empty/duplicate matrix identifiers in {path}; resolve upstream.")
    t = t.set_index(key)
    t.columns = t.columns.str.strip()
    if t.columns.duplicated().any():
        raise ValueError(f"Duplicate columns in {path}")
    t = t.apply(pd.to_numeric, errors="raise")
    if not np.isfinite(t.to_numpy()).all():
        raise ValueError(f"Missing/nonfinite matrix values in {path}; never zero-fill missing genes.")
    return t.astype(float)


def keyed_table(path, candidates, keys):
    """Select an explicit, unique identifier column matching every requested key.

    Supports both numbered official components and the named PRECISE411 exports
    used by the attached R script. Ambiguous/nonmatching files fail explicitly.
    """
    t = read_table(path)
    for col in list(dict.fromkeys(candidates + [t.columns[0]])):
        if col not in t:
            continue
        values = t[col].str.strip()
        if not values.duplicated().any() and set(keys).issubset(set(values)):
            t[col] = values
            return t.set_index(col).loc[keys]
    raise ValueError(f"No unique identifier column matches all {len(keys)} keys in {path}. "
                     "Use matching M/A/annotation files from the same reference release.")


def bh(p):
    p = np.asarray(p, float)
    out = np.full(p.shape, np.nan)
    valid = np.flatnonzero(np.isfinite(p))
    if len(valid):
        order = valid[np.argsort(p[valid])]
        q = p[order] * len(order) / np.arange(1, len(order) + 1)
        out[order] = np.minimum(1., np.minimum.accumulate(q[::-1])[::-1])
    return out


def save_table(t, path, index=True):
    t.to_csv(path, sep="\t", index=index, na_rep="NA")


def slug(s):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s))


def save_fig(fig, folder, stem):
    fig.savefig(folder / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(folder / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def label_points(ax, points, fontsize=6):
    """Greedy screen-space label placement for densely grouped PCA replicates."""
    ax.figure.canvas.draw()
    renderer = ax.figure.canvas.get_renderer()
    occupied = []
    candidates = [(6, 6), (6, -10), (-6, 6), (-6, -10)]
    for distance in [18, 30, 45, 60, 80, 105, 130]:
        candidates.extend([(6, distance), (6, -distance), (-6, distance), (-6, -distance),
                           (distance, 6), (distance, -10), (-distance, 6), (-distance, -10)])
    for name, x, y in points:
        artist = None
        for dx, dy in candidates:
            artist = ax.annotate(name, (x, y), fontsize=fontsize, xytext=(dx, dy), textcoords="offset points",
                                 ha="left" if dx >= 0 else "right", va="bottom" if dy >= 0 else "top")
            bbox = artist.get_window_extent(renderer).expanded(1.04, 1.12)
            if not any(bbox.overlaps(previous) for previous in occupied):
                occupied.append(bbox)
                break
            artist.remove()
            artist = None
        if artist is None:
            artist = ax.annotate(name, (x, y), fontsize=fontsize, xytext=(6, 150), textcoords="offset points")
        if abs(dx) + abs(dy) > 25:
            ax.annotate("", (x, y), xytext=(dx, dy), textcoords="offset points",
                        arrowprops={"arrowstyle": "-", "lw": .4, "color": "#AAAAAA"})


def projection(M, X, rcond):
    """SVD pseudoinverse with explicit numerical rank diagnostics."""
    u, s, vh = np.linalg.svd(M.to_numpy(), full_matrices=False)
    rank = int(np.sum(s > rcond * s[0]))
    if rank < M.shape[1]:
        raise ValueError(f"Shared-gene M has rank {rank}, fewer than {M.shape[1]} modules. "
                         "Improve gene overlap or check the reference.")
    pinv = (vh.T / s) @ u.T
    A = pd.DataFrame(pinv @ X.to_numpy(), index=M.columns, columns=X.columns)
    fitted = pd.DataFrame(M.to_numpy() @ A.to_numpy(), index=X.index, columns=X.columns)
    return A, fitted, s


def design_matrix(meta, cfg):
    levels = list(dict.fromkeys(meta["group"]))
    arrays = [(meta["group"] == g).to_numpy(float) for g in levels]
    names = ["group:" + g for g in levels]
    categorical = list(cfg.get("categorical_covariates", []))
    center = cfg.get("center_column")
    if center and center not in categorical:
        categorical.append(center)
    for col in categorical:
        if col == "group" or col not in meta or (meta[col].str.strip() == "").any():
            raise ValueError(f"Invalid or missing categorical covariate {col}")
        vals = list(dict.fromkeys(meta[col]))
        for v in vals[1:]:
            arrays.append((meta[col] == v).to_numpy(float))
            names.append(col + ":" + v)
    for col in cfg.get("numeric_covariates", []):
        if col not in meta:
            raise ValueError(f"Missing numeric covariate {col}")
        v = pd.to_numeric(meta[col], errors="raise").to_numpy(float)
        if not np.isfinite(v).all() or np.std(v) == 0:
            raise ValueError(f"Nonfinite or constant numeric covariate {col}")
        arrays.append(v - np.mean(v))
        names.append(col)
    D = np.column_stack(arrays)
    if np.linalg.matrix_rank(D) < D.shape[1]:
        raise ValueError("Design is confounded/rank deficient. A batch effect cannot be "
                         "estimated if condition and batch coincide completely.")
    return D, names, levels


def fit_contrasts(Y, D, names, meta, contrasts):
    """OLS using biological samples as observations; all modules fitted jointly.

    Assumes independent residuals and a common residual variance within each
    module. Fixed block terms can represent a matched design. A technical
    replicate must not be entered as an independent biological sample.
    """
    inv = np.linalg.inv(D.T @ D)
    beta = inv @ D.T @ Y.T.to_numpy()
    residual = Y.T.to_numpy() - D @ beta
    dof = D.shape[0] - D.shape[1]
    var = np.sum(residual ** 2, axis=0) / dof if dof > 0 else np.full(Y.shape[0], np.nan)
    rows = []
    for spec in contrasts:
        weights = spec["weights"]
        if not np.isclose(sum(weights.values()), 0):
            raise ValueError(f"Contrast {spec['id']} weights must sum to zero.")
        c = np.zeros(D.shape[1])
        ns = {}
        for group, weight in weights.items():
            if "group:" + group not in names:
                raise ValueError(f"Contrast {spec['id']} references absent group {group}")
            c[names.index("group:" + group)] = float(weight)
            ns[group] = int((meta.group == group).sum())
        if not any(c):
            raise ValueError("Contrast cannot be all zero.")
        effect = c @ beta
        se = np.sqrt(np.maximum(0, var * (c @ inv @ c)))
        testable = dof > 0 and all(n >= 2 for n in ns.values())
        p = np.full(len(effect), np.nan)
        good = (se > 1e-12) & np.isfinite(se)
        if testable:
            p[good] = 2 * stats.t.sf(np.abs(effect[good] / se[good]), dof)
        ci = stats.t.ppf(.975, dof) * se if testable else np.full(len(effect), np.nan)
        for j, module in enumerate(Y.index):
            rows.append({"contrast_id": spec["id"], "contrast_label": spec.get("label", spec["id"]),
                         "module_id": module, "delta_activity": effect[j], "standard_error": se[j],
                         "ci95_low": effect[j] - ci[j], "ci95_high": effect[j] + ci[j],
                         "p_value": p[j], "residual_df": dof,
                         "replicate_counts": json.dumps(ns, sort_keys=True),
                         "test_status": "tested" if testable and good[j] else "descriptive_only"})
    return pd.DataFrame(rows), residual, beta


def correlations(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    pearson = np.corrcoef(x, y)[0, 1] if np.std(x) > 0 and np.std(y) > 0 else np.nan
    spearman = stats.spearmanr(x, y).statistic if np.std(x) > 0 and np.std(y) > 0 else np.nan
    norm = np.linalg.norm(x) * np.linalg.norm(y)
    return pearson, spearman, (x @ y / norm if norm else np.nan)


def reference_analysis(cfg, base, M, ann, delta, out, min_effect, fdr):
    """Exploratory atlas matches + historical log-normal reference-null sensitivity.

    Uses every within-project/condition biological replicate pair. Unlike calling
    infer_activities on the new data, this procedure uses the atlas A for its null.
    It cannot validate the new experiment's technical variance. No reference
    samples are merged into the primary OLS design or primary PCA.
    """
    if not cfg.get("reference_A"):
        return
    RA = matrix(base / cfg["reference_A"], cfg.get("reference_A_index_column"))
    if not set(M.columns).issubset(RA.index):
        raise ValueError("Reference A rows do not match M columns.")
    RA = RA.loc[M.columns]
    rm = keyed_table(base / cfg["reference_samples"],
                     ["unique_identifier", "sample_id", "sample", "Run"], list(RA.columns))
    profiles = RA.copy()
    kind = "atlas_sample_vs_its_native_project_reference_EXPLORATORY"
    if cfg.get("reference_contrasts"):
        rc = read_table(base / cfg["reference_contrasts"])
        needed = {"contrast_id", "control_samples", "test_samples"}
        if not needed.issubset(rc):
            raise ValueError(f"Reference contrasts need {needed}")
        profiles = pd.DataFrame(index=RA.index)
        records = []
        for _, row in rc.iterrows():
            ctl = row.control_samples.split(";")
            tst = row.test_samples.split(";")
            if row.contrast_id in profiles:
                raise ValueError("Duplicate reference contrast IDs")
            if not set(ctl + tst).issubset(RA.columns) or set(ctl) & set(tst):
                raise ValueError(f"Invalid atlas samples in {row.contrast_id}")
            profiles[row.contrast_id] = RA[tst].mean(axis=1) - RA[ctl].mean(axis=1)
            records.append(row.to_dict())
        rm_profiles = pd.DataFrame(records).set_index("contrast_id")
        kind = "explicit_atlas_condition_minus_matched_control_EXPLORATORY"
    else:
        rm_profiles = rm
    matches = []
    for contrast in delta:
        for profile in profiles:
            p, s, c = correlations(delta[contrast], profiles[profile])
            matches.append({"contrast_id": contrast, "atlas_profile": profile,
                            "pearson_r": p, "spearman_rho": s, "cosine_similarity": c,
                            "comparison_type": kind, **rm_profiles.loc[profile].to_dict()})
    save_table(pd.DataFrame(matches).sort_values(["contrast_id", "pearson_r"], ascending=[True, False]),
               out / "atlas_profile_matches_EXPLORATORY.tsv", False)
    if not cfg.get("reference_null_dima", False):
        return
    project = cfg.get("reference_project_column", "project_name")
    condition = cfg.get("reference_condition_column", "condition")
    if project not in rm or condition not in rm:
        raise ValueError("Reference-null DiMA requires project and condition columns.")
    if (rm[[project, condition]].apply(lambda x: x.str.strip()) == "").any().any():
        raise ValueError("Reference project/condition fields contain blanks.")
    pairs = []
    for _, sub in rm.groupby([project, condition], sort=False):
        for s1, s2 in combinations(sub.index, 2):
            pairs.append(np.abs(RA[s1].to_numpy() - RA[s2].to_numpy()))
    if len(pairs) < 10:
        raise ValueError("Too few reference replicate pairs for module-specific log-normal nulls.")
    null = np.column_stack(pairs)
    save_table(pd.DataFrame(null, index=RA.index), out / "reference_within_condition_absolute_differences.tsv")
    fits, rows = [], []
    for j, module in enumerate(RA.index):
        values = null[j][null[j] > 0]
        params = None
        if len(values) >= 10 and np.std(np.log(values)) > 1e-12:
            # Explicit adaptation: two-parameter positive log-normal (loc=0).
            # The documented PyModulon util uses an unconstrained three-parameter fit.
            params = stats.lognorm.fit(values, floc=0)
        fits.append({"module_id": module, "n_pairs": len(pairs), "n_positive_pairs": len(values),
                     "shape": params[0] if params else np.nan, "loc": 0,
                     "scale": params[2] if params else np.nan,
                     "fit_status": "fitted_loc_zero" if params else "not_fit"})
        for contrast in delta:
            difference = delta.loc[module, contrast]
            p = stats.lognorm.sf(abs(difference), *params) if params else np.nan
            rows.append({"module_id": module, "contrast_id": contrast,
                         "delta_activity": difference, "p_reference_null": p})
    rr = pd.DataFrame(rows)
    rr["q_reference_null_per_contrast"] = rr.groupby("contrast_id")["p_reference_null"].transform(bh)
    rr["reference_null_DiMA"] = (rr.q_reference_null_per_contrast < fdr) & (rr.delta_activity.abs() > min_effect)
    rr = rr.join(ann, on="module_id", rsuffix="_annotation")
    save_table(rr, out / "DiMA_reference_null_SENSITIVITY.tsv", False)
    save_table(pd.DataFrame(fits), out / "reference_null_log_normal_fit.tsv", False)


def make_plots(A, results, delta, ann, meta, cfg, plots, gene_effects, M, membership):
    """PNG/PDF plots, with source numerical tables exported by the main pipeline."""
    labels = ann["display_name"].to_dict()
    contrast_labels = {c["id"]: c.get("label", c["id"]) for c in cfg["contrasts"]}
    sig = results.pivot(index="module_id", columns="contrast_id", values="DiMA").reindex_like(delta).fillna(False)
    selected = delta.index[sig.any(axis=1)]
    exploratory = len(selected) == 0
    if exploratory:
        selected = delta.abs().max(axis=1).nlargest(cfg.get("heatmap_max_modules", 40)).index
    else:
        selected = delta.loc[selected].abs().max(axis=1).nlargest(cfg.get("heatmap_max_modules", 40)).index
    heat = delta.loc[selected]
    if len(heat) > 1:
        heat = heat.iloc[leaves_list(linkage(heat.to_numpy(), method="average", metric="euclidean"))]
    vmax = max(float(heat.abs().to_numpy().max()), 1e-6)
    fig, ax = plt.subplots(figsize=(max(6, len(delta.columns) * 1.4), max(4, len(heat) * .25)))
    stars = sig.loc[heat.index, heat.columns].replace({True: "*", False: ""}).to_numpy()
    sns.heatmap(heat, cmap="RdBu_r", center=0, vmin=-vmax, vmax=vmax, ax=ax,
                annot=stars, fmt="", cbar_kws={"label": "Differential inferred activity (source units)"})
    ax.set_yticklabels([labels[m] for m in heat.index], rotation=0)
    ax.set_xticklabels([contrast_labels[c] for c in heat.columns], rotation=35, ha="right")
    ax.set(xlabel="", ylabel="", title=("Largest shifts (exploratory; no DiMAs)" if exploratory else
                                        "Differential iModulon activity; * primary DiMA"))
    save_fig(fig, plots, "delta_activity_heatmap")
    save_table(heat, plots.parent / "heatmap_delta_activity_source.tsv")
    # Model effects and uncertainty are more informative than a scatter whose
    # control mean is forced to zero by the reference-centering operation.
    ci_modules = list(heat.index[:min(16, len(heat))])
    fig, axes = plt.subplots(1, len(delta.columns), figsize=(max(6, len(delta.columns) * 2.5), max(4, len(ci_modules) * .35)),
                             sharey=True, squeeze=False)
    for ax, c in zip(axes[0], delta.columns):
        sub = results[results.contrast_id == c].set_index("module_id").loc[ci_modules]
        yy = np.arange(len(ci_modules))
        ax.axvline(0, color="#777777", lw=.8)
        ax.hlines(yy, sub.ci95_low, sub.ci95_high, color="#777777", lw=1)
        ax.scatter(sub.delta_activity, yy, c=np.where(sub.DiMA, "#B2182B", "#777777"), s=22)
        ax.set(title=contrast_labels[c], xlabel="Delta activity (95% CI)", yticks=yy)
        ax.set_yticklabels([labels[m] for m in ci_modules], fontsize=8)
    axes[0][0].invert_yaxis()
    fig.tight_layout()
    save_fig(fig, plots, "delta_activity_model_confidence_intervals")
    # A control-centered mean A is approximately zero. These scatter plots are
    # conventional displays of component differences, not independent evidence.
    groups = meta.groupby("group", sort=False).groups
    means = pd.DataFrame({g: A[list(samples)].mean(axis=1) for g, samples in groups.items()})
    simple = [c for c in cfg["contrasts"] if len(c["weights"]) == 2 and sorted(c["weights"].values()) == [-1, 1]]
    for c in simple:
        ctl = next(g for g, w in c["weights"].items() if w == -1)
        tst = next(g for g, w in c["weights"].items() if w == 1)
        hit = sig[c["id"]]
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.scatter(means[ctl], means[tst], c=np.where(hit, "#B2182B", "#9B9B9B"), s=22)
        lim = np.max(np.abs(means[[ctl, tst]].to_numpy())) * 1.15 + 1
        xx = np.array([-lim, lim])
        for shift, style in [(0, "-"), (cfg.get("min_abs_delta", 5), "--"), (-cfg.get("min_abs_delta", 5), "--")]:
            ax.plot(xx, xx + shift, style, color="#777777", lw=.8)
        for m in delta.index[hit][:8]:
            ax.annotate(labels[m], (means.loc[m, ctl], means.loc[m, tst]), fontsize=7)
        ax.set(xlabel=f"{ctl}: mean inferred activity", ylabel=f"{tst}: mean inferred activity",
               title=c.get("label", c["id"]), xlim=(-lim, lim), ylim=(-lim, lim))
        save_fig(fig, plots, "DiMA_scatter_" + slug(c["id"]))
    # Select a few modules for replicate plots + conditional gene loading checks.
    priority = results.sort_values(["DiMA", "q_global", "delta_activity"], ascending=[False, True, False])
    top = list(dict.fromkeys(priority.module_id))[:cfg.get("validation_top_modules", 6)]
    custom = cfg.get("validation_modules", [])
    for m in custom:
        if m not in A.index:
            raise ValueError(f"Requested validation module {m} absent from M")
    top = list(dict.fromkeys(custom + top))
    counts = meta.group.value_counts()
    for m in top:
        tab = meta.copy()
        tab["activity"] = A.loc[m]
        fig, ax = plt.subplots(figsize=(max(7, counts.size * 1.1), 4))
        sns.boxplot(data=tab, x="group", y="activity", order=list(groups), color="#E5E5E5", ax=ax, fliersize=0)
        sns.stripplot(data=tab, x="group", y="activity", order=list(groups), hue="batch" if "batch" in tab else None,
                      dodge=False, jitter=.08, ax=ax, size=5)
        ax.set(title=labels[m], xlabel="", ylabel="Inferred activity (source orientation)")
        ax.tick_params(axis="x", rotation=30)
        save_fig(fig, plots, "replicate_activity_" + slug(m))
        for c in cfg["contrasts"]:
            # Remove contributions of all OTHER modules to inspect this component.
            total = M.to_numpy() @ delta[c["id"]].to_numpy()
            marginal = gene_effects[c["id"]].to_numpy()
            contribution = M[m].to_numpy() * delta.loc[m, c["id"]]
            conditional = marginal - (total - contribution)
            member = membership[m].to_numpy()
            fig, axes = plt.subplots(1, 2, figsize=(10, 4))
            for ax, yy, title in [(axes[0], marginal, "Observed log-expression contrast"),
                                  (axes[1], conditional, "Other fitted modules removed (descriptive)")]:
                ax.scatter(M[m], yy, c=np.where(member, "#B2182B", "#C5C5C5"), s=9, alpha=.6)
                xx = np.array([M[m].min(), M[m].max()])
                ax.plot(xx, xx * delta.loc[m, c["id"]], color="#2166AC", lw=1)
                ax.set(xlabel="Source M loading", ylabel="log2 expression difference", title=title)
            fig.suptitle(f"{labels[m]} | {contrast_labels[c['id']]}", fontsize=11)
            fig.tight_layout()
            save_fig(fig, plots, "loading_check_" + slug(m) + "_" + slug(c["id"]))
    # Unscaled activity PCA. Components retain the reference's M scaling.
    Q = A.T.to_numpy()
    Q = Q - Q.mean(axis=0)
    u, s, vh = np.linalg.svd(Q, full_matrices=False)
    if len(s) >= 2 and np.sum(s ** 2) > 0:
        scores = u[:, :2] * s[:2]
        pcs = pd.DataFrame(scores, index=A.columns, columns=["PC1", "PC2"]).join(meta)
        save_table(pcs, plots.parent / "activity_PCA_scores.tsv")
        save_table(pd.DataFrame(vh[:2].T, index=A.index, columns=["PC1", "PC2"]), plots.parent / "activity_PCA_loadings.tsv")
        pct = s ** 2 / np.sum(s ** 2) * 100
        fig, ax = plt.subplots(figsize=(8, 6))
        palette = dict(zip(groups, sns.color_palette("tab10", len(groups))))
        for g, sub in pcs.groupby("group", sort=False):
            ax.scatter(sub.PC1, sub.PC2, color=palette[g], label=g, s=35)
        centers = pcs.groupby("group")[["PC1", "PC2"]].mean()
        control = cfg["control_group"]
        for g in groups:
            if g != control:
                ax.annotate("", xy=centers.loc[g], xytext=centers.loc[control],
                            arrowprops={"arrowstyle": "->", "color": palette[g], "lw":1})
        ax.set(xlabel=f"PC1 ({pct[0]:.1f}%)", ylabel=f"PC2 ({pct[1]:.1f}%)",
               title="Inferred activity space (descriptive; unscaled)")
        ax.margins(.17)
        ax.legend(fontsize=8, bbox_to_anchor=(1.02, 1), loc="upper left")
        label_points(ax, [(sample, row.PC1, row.PC2) for sample, row in pcs.iterrows()])
        save_fig(fig, plots, "activity_PCA")


def run(config_path):
    config_path = Path(config_path).resolve()
    cfg = json.loads(config_path.read_text())
    base = config_path.parent
    out = (base / cfg["output_dir"]).resolve()
    out.mkdir(parents=True, exist_ok=True)
    plots = out / "figures"
    plots.mkdir(exist_ok=True)
    handler = logging.FileHandler(out / "run.log", mode="w")
    LOG.addHandler(handler)
    LOG.setLevel(logging.INFO)
    LOG.info("Pipeline version %s", VERSION)
    if cfg.get("input_units") != "log2_TPM":
        raise ValueError("Declare input_units='log2_TPM'. Provide sample-level log2(TPM+pseudocount), "
                         "not DEG tables, z scores, logCPM or raw counts.")
    if "input_pseudocount" not in cfg:
        raise ValueError("Declare input_pseudocount to document the upstream transformation.")
    pc = float(cfg["input_pseudocount"])
    if pc < 0:
        raise ValueError("Pseudocount must be nonnegative.")
    # 1. Load the expression matrix, fixed atlas and biological sample groups.
    # Use the broad gene-level TPM export, not a filtered DEG or logCPM table.
    X = matrix(base / cfg["expression"], cfg.get("expression_index_column"))
    if pc >= 1 and (X.to_numpy() < -1e-9).any():
        raise ValueError("Negative values contradict log2(TPM+pseudocount) with pseudocount >=1.")
    Mfull = matrix(base / cfg["M"], cfg.get("M_index_column"))
    if Mfull.shape[0] <= Mfull.shape[1]:
        raise ValueError("M must have genes as rows and iModulons as columns.")
    ann = keyed_table(base / cfg["imodulon_table"], ["imodulon_name", "name", "k", "module_id"], list(Mfull.columns))
    for col in ["display_name", "regulator", "category", "function_description"]:
        if col not in ann:
            ann[col] = ""
    namecol = next((c for c in ["imodulon_name", "name"] if c in ann), None)
    ann["display_name"] = ann[namecol].where(ann[namecol] != "", ann.index.to_series()) if namecol else ann.index
    ann.index.name = "module_id"
    meta = keyed_table(base / cfg["samples"], ["sample_id", "sample", "unique_identifier"], list(X.columns))
    if "group" not in meta or (meta.group.str.strip() == "").any():
        raise ValueError("Sample metadata needs a nonempty 'group' column.")
    # All X samples are analyzed. Remove technical replicates/excluded samples upstream.
    meta.index.name = "sample_id"
    control = cfg["control_group"]
    if control not in set(meta.group):
        raise ValueError(f"Control group {control} absent in metadata")
    contrasts = cfg["contrasts"]
    ids = [c["id"] for c in contrasts]
    if len(set(ids)) != len(ids) or not ids:
        raise ValueError("Contrasts must have unique nonempty IDs.")
    min_effect = float(cfg.get("min_abs_delta", 5))
    fdr = float(cfg.get("fdr", .05))
    if min_effect < 0 or not 0 < fdr < 1:
        raise ValueError("Invalid effect/FDR settings")
    # 2. Align expression and atlas by PAO1 locus tag, then center on controls.
    # Centering defines the reference state; expression is not z-scored.
    matched = [g for g in Mfull.index if g in X.index]
    missing = Mfull.index.difference(X.index)
    save_table(pd.DataFrame({"gene_id": missing}), out / "missing_reference_genes.tsv", False)
    save_table(pd.DataFrame({"gene_id": X.index.difference(Mfull.index)}), out / "expression_genes_not_in_M.tsv", False)
    fraction = len(matched) / len(Mfull)
    if fraction < cfg.get("min_gene_overlap_fraction", .9):
        raise ValueError(f"Only {fraction:.1%} reference genes matched. Check PAO1 locus tags.")
    M = Mfull.loc[matched]
    X = X.loc[matched]
    center = cfg.get("center_column")
    centered = X.copy()
    if center:
        if center not in meta or (meta[center].str.strip() == "").any():
            raise ValueError("Invalid center_column in sample metadata")
        for block, sub in meta.groupby(center, sort=False):
            controls = sub.index[sub.group == control]
            if len(controls) == 0:
                raise ValueError(f"No matched {control} controls in {center}={block}")
            centered.loc[:, sub.index] = X[sub.index].sub(X[controls].mean(axis=1), axis=0)
    else:
        centered = X.sub(X[meta.index[meta.group == control]].mean(axis=1), axis=0)
    # Preserve every source M column and its native sign/scale. No thresholds in
    # projection, no row z scoring, no sequential single-module projections.
    # 3. Infer every module jointly using the full signed loading matrix.
    # Membership thresholds describe module genes; they do not truncate M.
    A, reconstructed, singular = projection(M, centered, cfg.get("pinv_rcond", 1e-12))
    membership = pd.DataFrame(False, index=M.index, columns=M.columns)
    moduleqc = []
    thresholds = None
    if cfg.get("M_thresholds"):
        thresholds = keyed_table(base / cfg["M_thresholds"], ["imodulon_name", "module_id", "k", "imodulon_k"], list(M.columns))
        thresholds = pd.to_numeric(thresholds["threshold"], errors="raise")
        if not np.isfinite(thresholds).all() or (thresholds < 0).any():
            raise ValueError("Invalid membership thresholds")
    for module in M.columns:
        wfull = Mfull[module]
        w = M[module]
        retained = float(np.sum(w ** 2) / np.sum(wfull ** 2)) if np.sum(wfull ** 2) else 0
        top = w.abs().idxmax()
        members = wfull.abs() > thresholds[module] if thresholds is not None else None
        if members is not None:
            membership[module] = members.loc[M.index]
        moduleqc.append({"module_id": module, "retained_squared_loading_fraction": retained,
                         "top_loading_gene": top, "top_gene_squared_loading_fraction": float(w[top] ** 2 / np.sum(w ** 2)),
                         "n_reference_members": int(members.sum()) if members is not None else np.nan,
                         "n_shared_members": int(membership[module].sum()) if members is not None else np.nan,
                         "projection_coverage_ok": retained >= cfg.get("min_module_weight_fraction", .9)})
    mq = pd.DataFrame(moduleqc).set_index("module_id")
    # 4. Fit biological samples, rather than pre-averaged group activities.
    # With the supplied config: 18 samples, six groups and five contrasts.
    D, names, levels = design_matrix(meta, cfg)
    results, residual, beta = fit_contrasts(A, D, names, meta, contrasts)
    gene_res, _, _ = fit_contrasts(centered, D, names, meta, contrasts)
    gene_effects = gene_res.pivot(index="module_id", columns="contrast_id", values="delta_activity").reindex(index=M.index, columns=ids)
    # The manuscript uses GLOBAL BH over all 116 x 5 tests.
    # The per-contrast q column is retained for reference, not primary calls.
    results["q_per_contrast"] = results.groupby("contrast_id")["p_value"].transform(bh)
    results["q_global"] = bh(results.p_value)
    scope = cfg.get("fdr_scope", "global")
    if scope not in ["global", "per_contrast"]:
        raise ValueError("fdr_scope must be global or per_contrast")
    results["q_primary"] = results["q_" + scope]
    results = results.join(mq, on="module_id").join(ann, on="module_id", rsuffix="_annotation")
    results["large_effect"] = results.delta_activity.abs() > min_effect
    results["DiMA"] = (results.q_primary < fdr) & results.large_effect & results.projection_coverage_ok
    results["direction"] = np.where(results.delta_activity > 0, "higher_module_signature", "lower_module_signature")
    delta = results.pivot(index="module_id", columns="contrast_id", values="delta_activity").reindex(index=M.columns, columns=ids)
    # 5. Keep the complete table: downstream scripts need nonsignificant
    # modules too, not just the 29 calls in this manuscript dataset.
    save_table(A, out / "A_inferred_sample_activities.tsv")
    save_table(delta, out / "delta_A_model_contrasts.tsv")
    save_table(centered, out / "centered_log2_TPM_shared_genes.tsv")
    save_table(pd.DataFrame(D, index=meta.index, columns=names), out / "model_design_matrix.tsv")
    save_table(meta, out / "sample_metadata_used.tsv")
    save_table(results, out / "DiMA_primary_all.tsv", False)
    save_table(results[results.DiMA].sort_values("q_primary"), out / "DiMA_primary_significant.tsv", False)
    save_table(mq.join(ann, rsuffix="_annotation"), out / "module_projection_QC.tsv")
    # Do not fabricate requested regulators that the atlas did not recover.
    searchable = ann.astype(str).agg(" ".join, axis=1) + " " + ann.index.to_series()
    terms = cfg.get("stress_annotation_terms", ["AmgR", "AmgS", "CprR", "CprS", "OxyR", "SoxR", "Anr", "AlgU", "RpoH", "MexZ"])
    hits = []
    for term in terms:
        found = searchable.str.contains(re.escape(term), case=False, regex=True)
        matching = list(ann.index[found])
        if matching:
            hits.extend({"search_term": term, "module_id": m, "display_name": ann.loc[m, "display_name"],
                         "annotation_match": True} for m in matching)
        else:
            hits.append({"search_term": term, "module_id": "", "display_name": "",
                         "annotation_match": False})
    save_table(pd.DataFrame(hits), out / "candidate_stress_regulator_annotation_matches.tsv", False)
    save_table(gene_effects, out / "gene_log_expression_model_contrasts.tsv")
    save_table(pd.DataFrame({"singular_value": singular}), out / "M_shared_singular_values.tsv", False)
    # Projection goodness of fit uses centered sum of squares, not the original
    # atlas's explained variance. It can be negative for arbitrary reconstructions;
    # least-squares projection should be >=0 within numerical precision.
    rr = centered - reconstructed
    denom = np.sum(centered.to_numpy() ** 2, axis=0)
    sse = np.sum(rr.to_numpy() ** 2, axis=0)
    quality = pd.DataFrame({"projection_R2_centered": np.divide(denom - sse, denom, out=np.full_like(denom, np.nan), where=denom > 1e-12),
                            "residual_RMSE": np.sqrt(np.mean(rr.to_numpy() ** 2, axis=0))}, index=A.columns)
    save_table(quality.join(meta), out / "sample_projection_QC.tsv")
    save_table(pd.DataFrame(residual.T, index=A.index, columns=A.columns), out / "OLS_activity_residuals.tsv")
    save_table(centered.corr(), out / "centered_expression_sample_correlations.tsv")
    save_table(X.corr(), out / "log2TPM_sample_correlations.tsv")
    contrastqc, residualgene, validation = [], [], []
    for c in ids:
        observed = gene_effects[c].to_numpy()
        fit = M.to_numpy() @ delta[c].to_numpy()
        denom = observed @ observed
        contrastqc.append({"contrast_id": c, "projection_R2_centered": 1 - np.sum((observed - fit) ** 2) / denom if denom > 1e-12 else np.nan})
        order = np.argsort(np.abs(observed - fit))[::-1][:100]
        for i in order:
            residualgene.append({"contrast_id": c, "gene_id": M.index[i], "observed_log_expression_delta": observed[i],
                                 "reconstructed_delta": fit[i], "residual": observed[i] - fit[i]})
        for m in M.columns:
            member = membership[m].to_numpy()
            conditional = observed - (fit - M[m].to_numpy() * delta.loc[m, c])
            if member.sum() >= 3:
                r, _, _ = correlations(M[m].to_numpy()[member], observed[member])
                rc, _, _ = correlations(M[m].to_numpy()[member], conditional[member])
                concord = np.mean(np.sign(observed[member]) == np.sign(M[m].to_numpy()[member] * delta.loc[m, c]))
            else:
                r = rc = concord = np.nan
            validation.append({"module_id": m, "contrast_id": c, "member_loading_vs_observed_pearson": r,
                               "member_loading_vs_conditional_pearson_DESCRIPTIVE": rc,
                               "member_sign_concordance": concord})
    save_table(pd.DataFrame(contrastqc), out / "contrast_projection_QC.tsv", False)
    save_table(pd.DataFrame(residualgene), out / "largest_unexplained_gene_changes.tsv", False)
    save_table(pd.DataFrame(validation), out / "module_gene_coherence_DESCRIPTIVE.tsv", False)
    # Shared vs specific: significance is a detection label; nonsignificance in
    # another condition does not prove specificity. Interaction contrasts needed.
    calls = results.pivot(index="module_id", columns="contrast_id", values="DiMA").reindex(index=M.columns, columns=ids)
    calls["n_detected_contrasts"] = calls.sum(axis=1)
    calls["detected_contrasts"] = calls[ids].apply(lambda r: ";".join(r.index[r]), axis=1)
    save_table(calls.join(ann), out / "shared_detected_module_signatures.tsv")
    tob = cfg.get("tobramycin_contrast")
    if tob:
        if tob not in delta:
            raise ValueError("tobramycin_contrast absent from configured contrasts")
        compare = []
        for c in ids:
            if c == tob:
                continue
            p, s, cos = correlations(delta[c], delta[tob])
            both = calls[c] & calls[tob]
            same = np.sign(delta[c]) == np.sign(delta[tob])
            compare.append({"contrast_id": c, "tobramycin_contrast": tob, "pearson_r": p,
                            "spearman_rho": s, "cosine_similarity": cos,
                            "n_both_DiMA_same_direction": int((both & same).sum()),
                            "n_both_DiMA_opposite_direction": int((both & ~same).sum()),
                            "interpretation": "descriptive_shared_state_not_causal_sensitization"})
        save_table(pd.DataFrame(compare), out / "knockdown_vs_tobramycin_similarity_DESCRIPTIVE.tsv", False)
        fig, ax = plt.subplots(figsize=(6, 4))
        ct = pd.DataFrame(compare)
        ax.bar(ct.contrast_id, ct.cosine_similarity, color="#2166AC")
        ax.axhline(0, color="black", lw=.7)
        ax.set(ylim=(-1, 1), ylabel="Cosine similarity to TOB activity change", title="Shared activity pattern (descriptive)")
        ax.tick_params(axis="x", rotation=35)
        save_fig(fig, plots, "knockdown_TOB_profile_similarity")
    # Optional sensitivity refit: excluding the direct CRISPRi target genes
    # helps assess whether activity shifts are dominated by these genes.
    excludes = cfg.get("sensitivity_exclude_gene_ids", [])
    if excludes:
        unknown = set(excludes) - set(M.index)
        if unknown:
            raise ValueError(f"Sensitivity exclusions not among shared genes: {sorted(unknown)}")
        genes = M.index.difference(excludes, sort=False)
        As, _, _ = projection(M.loc[genes], centered.loc[genes], cfg.get("pinv_rcond", 1e-12))
        rs, _, _ = fit_contrasts(As, D, names, meta, contrasts)
        rs["q_global_excluded"] = bh(rs.p_value)
        rs = rs.rename(columns={"delta_activity": "delta_activity_excluded", "p_value": "p_excluded"})
        rs = rs[["module_id", "contrast_id", "delta_activity_excluded", "p_excluded", "q_global_excluded"]]
        rs = results[["module_id", "contrast_id", "delta_activity", "DiMA"]].merge(rs, on=["module_id", "contrast_id"])
        rs["same_effect_direction"] = np.sign(rs.delta_activity) == np.sign(rs.delta_activity_excluded)
        rs["effect_change"] = rs.delta_activity_excluded - rs.delta_activity
        save_table(rs, out / "excluded_target_gene_projection_SENSITIVITY.tsv", False)
    reference_analysis(cfg, base, M, ann, delta, out, min_effect, fdr)
    if cfg.get("plots", True):
        make_plots(A, results, delta, ann, meta, cfg, plots, gene_effects, M, membership)
    messages = ["Activities are inferred transcriptomic signatures, not direct measurements of TF activity.",
                "Primary test assumes independent biological replicates and OLS residual assumptions.",
                "No significance is assigned if any contributing group has <2 replicates or residual variance is degenerate.",
                "TOB similarity uses a shared control; correlation is descriptive and is not a synergy test.",
                "Negative/sign-switched ICA components require gene-level interpretation; a repressor-associated module may increase when its regulator is less active."]
    if thresholds is None:
        messages.append("No membership thresholds supplied: member-based coherence metrics unavailable.")
    if np.max(singular) / np.min(singular) > 1e4:
        messages.append("M is poorly conditioned; activities may be numerically unstable.")
    if (~mq.projection_coverage_ok).any():
        messages.append("Poorly covered modules excluded from primary DiMA calls; inspect module_projection_QC.tsv.")
    for message in messages:
        LOG.warning(message)
    (out / "interpretation_notes.txt").write_text("\n".join(messages) + "\n")
    # 6. Record input hashes and installed versions so a run can be traced.
    files = ["expression", "samples", "M", "imodulon_table", "M_thresholds", "reference_A", "reference_samples", "reference_contrasts"]
    hashes = {k: hashlib.sha256((base / cfg[k]).read_bytes()).hexdigest() for k in files if cfg.get(k)}
    provenance = {"pipeline_version": VERSION, "config": cfg, "input_sha256": hashes,
                  "matched_genes": len(matched), "reference_genes": len(Mfull), "gene_overlap_fraction": fraction,
                  "n_modules": len(M.columns), "M_rank": len(singular), "M_condition_number": float(singular[0] / singular[-1]),
                  "statistics": "sample-level OLS, two-sided t, BH; reference null is separate",
                  "versions": {"python": sys.version, "numpy": np.__version__, "pandas": pd.__version__,
                               "scipy": __import__("scipy").__version__, "matplotlib": matplotlib.__version__, "seaborn": sns.__version__}}
    (out / "run_provenance.json").write_text(json.dumps(provenance, indent=2))
    LOG.info("Finished: %d samples, %d shared genes, %d modules, %d primary DiMAs. Outputs: %s",
             X.shape[1], len(M), M.shape[1], int(results.DiMA.sum()), out)
    LOG.removeHandler(handler)
    handler.close()
    # Make the hand-off to the R figure/table scripts visible at completion.
    print(f"Underlying iModulon results: {out / 'DiMA_primary_all.tsv'}")
    print(f"Sample activities: {out / 'A_inferred_sample_activities.tsv'}")
    print(f"Projection checks: {out / 'module_projection_QC.tsv'}")
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", required=True, help="JSON configuration; relative paths resolved beside this file")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    try:
        run(args.config)
    except (ValueError, KeyError, FileNotFoundError, np.linalg.LinAlgError) as exc:
        LOG.error("Stopped: %s", exc)
        raise SystemExit(2)


if __name__ == "__main__":
    main()
