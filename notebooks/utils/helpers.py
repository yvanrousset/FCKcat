import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde, pearsonr
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from sklearn.metrics import r2_score, matthews_corrcoef

__all__ = [
    'compute_delta_kcat',
    'unpack_pairs_df',
    'compute_average_kcat',
    'plot_variance_decomposition',
    'compute_average_delta',
    'order_pairs_mutant_first',
    '_to_numeric_list',
    '_prep_xy',
    '_scatter_panel',
    '_scatter_panel_v2',
    '_shared_row_limits',
    '_metrics_no_threshold',
    '_bars_panel',
    'add_esm2',
    'add_esm2_max',
    'add_rxnfp',
    'add_symetric_pairs',
    'prepro',
    'prepro_v2',
    'add_esm1v',
]


def compute_delta_kcat(df):
    """Add a delta column = kcat - wildtype kcat, per Identifier group. Returns a copy."""
    df = df.copy()
    df['delta'] = None
    for ind in df.index:
        identifier = df.at[ind, "Identifier"]
        df_temp = df[df['Identifier'] == identifier]
        df_temp_wt = df_temp[df_temp['Type'] == 'wildtype']
        try:
            kcat_wt = df_temp_wt['geom_kcat_log10'].values[0]
            df.at[ind, 'delta'] = df.at[ind, 'geom_kcat_log10'] - kcat_wt
        except:
            df.at[ind, 'delta'] = None
    df["delta"] = pd.to_numeric(df["delta"], errors="coerce")
    return df


def unpack_pairs_df(df_pairs):
    """Unpack a paired DataFrame (list-valued columns) into individual rows, de-duplicated."""
    left_rows, right_rows = [], []
    for ind in df_pairs.index:
        base = {'Identifier': df_pairs.loc[ind, 'Identifier']}
        left_rows.append({**base,
                          'Type':      df_pairs.loc[ind, 'Type'][0],
                          'index':     df_pairs.loc[ind, 'index'][0],
                          'log10kcat': df_pairs.loc[ind, 'log10kcat'][0]})
        right_rows.append({**base,
                           'Type':      df_pairs.loc[ind, 'Type'][1],
                           'index':     df_pairs.loc[ind, 'index'][1],
                           'log10kcat': df_pairs.loc[ind, 'log10kcat'][1]})
    return (
        pd.concat([pd.DataFrame(left_rows), pd.DataFrame(right_rows)], ignore_index=True)
          .drop_duplicates()
          .reset_index(drop=True)
    )


def _to_numeric_list(val):
    """Return a flat list of numeric values from val (handles scalars and sequences)."""
    if isinstance(val, (list, tuple, np.ndarray, pd.Series)):
        return pd.to_numeric(pd.Series(val), errors='coerce').dropna().tolist()
    return pd.to_numeric(pd.Series([val]), errors='coerce').dropna().tolist()


def compute_average_kcat(df_test_unpaired, df_train_unpaired):
    """Compute average_kcat_mut and average_kcat_group for each test row from training data."""
    df = df_test_unpaired.copy()
    df['average_kcat_mut']   = np.nan
    df['average_kcat_group'] = np.nan

    for ind in df.index:
        groupid = df.loc[ind, 'Identifier']
        df_temp = df_train_unpaired[df_train_unpaired['Identifier'] == groupid]

        if df_temp.empty:
            df_temp_test_wt = df[(df['Identifier'] == groupid) & (df['Type'] == 'wildtype')]
            if not df_temp_test_wt.empty:
                df.at[ind, 'average_kcat_group'] = df_temp_test_wt['log10kcat'].values[0]
            continue

        mutant_vals, all_vals = [], []
        for i in df_temp.index:
            vals_i = _to_numeric_list(df_temp.loc[i, 'log10kcat'])
            all_vals.extend(vals_i)
            if df_temp.loc[i, 'Type'] == 'mutant':
                mutant_vals.extend(vals_i)

        df.at[ind, 'average_kcat_mut']   = np.mean(mutant_vals) if mutant_vals else np.nan
        df.at[ind, 'average_kcat_group'] = np.mean(all_vals)    if all_vals    else np.nan

    return df


def plot_variance_decomposition(between_group, within_group, labels, ax=None, fontsize=12):
    """Stacked bar chart of between/within group variance fractions with ICC annotations."""
    high_contrast = ['#004488', '#DDAA33', '#BB5566', '#000000']
    total         = between_group + within_group
    between_frac  = between_group / total
    within_frac   = within_group  / total
    icc           = between_frac

    if ax is None:
        fig, ax = plt.subplots(figsize=(4.2, 4.8))
    else:
        fig = ax.get_figure()

    x = np.arange(len(labels))
    width = 0.72
    ax.bar(x, between_frac, width=width, color=high_contrast[0], label='Between-group')
    ax.bar(x, within_frac,  width=width, bottom=between_frac,    color=high_contrast[1], label='Within-group')

    for i in range(len(labels)):
        ax.text(x[i], between_frac[i] / 2,
                f'{100*between_frac[i]:.1f}%', ha='center', va='center', color='white',  fontsize=fontsize)
        ax.text(x[i], between_frac[i] + within_frac[i] / 2,
                f'{100*within_frac[i]:.1f}%',  ha='center', va='center', color='black', fontsize=fontsize)
        ax.text(x[i], 1.02, f'ICC = {icc[i]:.2f}', ha='center', va='bottom',            fontsize=fontsize)

    ax.set_ylabel('Fraction of variance', fontsize=fontsize + 2)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=fontsize)
    ax.set_ylim(0, 1.12)
    ax.set_yticks(np.linspace(0, 1, 6))
    ax.legend(frameon=False, loc='lower center', bbox_to_anchor=(0.5, 1.08), ncol=2,
              fontsize=fontsize)
    return fig, ax


def compute_paired_delta_from_training(df_test, df_test_unpaired):
    """Add average_delta and average_delta_with_zero to df_test, computed from training averages."""
    df = df_test.copy()
    df['average_delta']           = np.nan
    df['average_delta_with_zero'] = np.nan

    for ind in df.index:
        types   = df.loc[ind, 'Type']
        groupid = df.loc[ind, 'Identifier']
        df_temp = df_test_unpaired[df_test_unpaired['Identifier'] == groupid]

        if df_temp.empty:
            df.at[ind, 'average_delta']           = 0
            df.at[ind, 'average_delta_with_zero'] = np.nan
        else:
            kcat_mutant = df_temp[df_temp['Type'] == 'mutant']['average_kcat_mut'].to_list()[0]
            if types[0] == 'wildtype' and types[1] == 'mutant':
                kcat_wt = df.loc[ind, 'log10kcat'][0]
                df.at[ind, 'average_delta']           = kcat_wt - kcat_mutant
                df.at[ind, 'average_delta_with_zero'] = kcat_wt - kcat_mutant
            elif types[1] == 'wildtype' and types[0] == 'mutant':
                kcat_wt = df.loc[ind, 'log10kcat'][1]
                df.at[ind, 'average_delta']           = kcat_mutant - kcat_wt
                df.at[ind, 'average_delta_with_zero'] = kcat_mutant - kcat_wt
    return df


def compute_group_average_delta(df):
    """
    Add an average_delta column per Identifier group.

    average_delta = mean(geom_kcat_log10 of mutants in group)
                    - geom_kcat_log10 of wildtype in group

    Returns
    -------
    df_out : pd.DataFrame
        Copy of input dataframe with a new 'average_delta' column.
    avg_delta_dict : dict
        Dictionary {Identifier: average_delta}, one value per group.
    """
    df_out = df.copy()
    df_out["average_delta"] = np.nan

    avg_delta_dict = {}

    for identifier, df_group in df_out.groupby("Identifier", sort=False):
        wt_vals = pd.to_numeric(
            df_group.loc[df_group["Type"] == "wildtype", "geom_kcat_log10"],
            errors="coerce"
        ).dropna()

        mut_vals = pd.to_numeric(
            df_group.loc[df_group["Type"] == "mutant", "geom_kcat_log10"],
            errors="coerce"
        ).dropna()

        if len(wt_vals) == 0 or len(mut_vals) == 0:
            avg_delta = np.nan
        else:
            # if there is exactly one wildtype per group, this is fine
            # if there can be several, you may want wt_vals.mean() instead
            wt_val = wt_vals.iloc[0]
            avg_delta = mut_vals.mean() - wt_val

        df_out.loc[df_group.index, "average_delta"] = avg_delta
        avg_delta_dict[identifier] = avg_delta

    return df_out, avg_delta_dict

def order_pairs_mutant_first(df):
    """Return a copy of df with every pair reordered so mutant is index 0, wildtype index 1."""
    df = df.copy()
    for ind in df.index:
        types = df.loc[ind, 'Type']
        if types[0] == 'mutant' and types[1] == 'wildtype':
            continue
        elif types[0] == 'wildtype' and types[1] == 'mutant':
            df.at[ind, 'Type']                    = [types[1], types[0]]
            df.at[ind, 'index']                   = [df.at[ind, 'index'][1],    df.at[ind, 'index'][0]]
            df.at[ind, 'log10kcat']               = [df.at[ind, 'log10kcat'][1], df.at[ind, 'log10kcat'][0]]
            df.at[ind, 'average_delta']           = -df.at[ind, 'average_delta']
            df.at[ind, 'average_delta_with_zero'] = -df.at[ind, 'average_delta_with_zero']
    return df


def _prep_xy(x_series, y_series):
    """Align two series, cast to float, drop non-finite values."""
    s_x = pd.to_numeric(x_series, errors="coerce")
    s_y = pd.to_numeric(y_series, errors="coerce")
    s_x, s_y = s_x.align(s_y, join="inner")
    x = s_x.to_numpy(dtype=float)
    y = s_y.to_numpy(dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    return x[mask], y[mask]


def _scatter_panel(ax, x_series, y_series, title, x_label, y_label,
                   label_fs=18, title_fs=16, annot_fs=20, cbar_fs=14):
    """Density scatter plot with 1:1 line, R² annotation, and inset colorbar."""
    x, y = _prep_xy(x_series, y_series)
    if x.size == 0:
        return None, None

    r2  = r2_score(x, y) if x.size >= 2 else np.nan
    xy  = np.vstack([x, y])
    try:
        z = gaussian_kde(xy)(xy)
    except Exception:
        z = np.ones_like(x)

    idx    = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]
    norm   = mpl.colors.Normalize(vmin=float(np.min(z)), vmax=float(np.max(z)))

    sc = ax.scatter(x, y, c=z, s=70, cmap="viridis", norm=norm, alpha=0.8, edgecolor="none")

    lo, hi = float(min(x.min(), y.min())), float(max(x.max(), y.max()))
    ax.plot([lo, hi], [lo, hi], "--", color="black", linewidth=1.0)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel(x_label, fontsize=label_fs, labelpad=10)
    ax.set_ylabel(y_label, fontsize=label_fs, labelpad=10)
    ax.set_title(title, fontsize=title_fs, weight="bold")
    ax.text(
        0.05, 0.95,
        rf"$R^2 = {r2:.2f}$" if np.isfinite(r2) else r"$R^2$ n/a",
        transform=ax.transAxes, fontsize=annot_fs, fontweight="bold", va="top",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                  edgecolor="black", linewidth=1.2, alpha=0.9)
    )

    axins = inset_axes(ax, width="35%", height="5%", loc="lower right", borderpad=2.0)
    cbar  = ax.figure.colorbar(sc, cax=axins, orientation="horizontal")
    cbar.ax.set_in_layout(False)
    cbar.set_label("Density", fontsize=cbar_fs, labelpad=10)
    cbar.ax.xaxis.set_label_position("top")
    cbar.ax.tick_params(labelsize=cbar_fs - 2)

    return x, y


def _scatter_panel_v2(ax, x_series, y_series, title, x_label, y_label):
    """Like _scatter_panel but renders a 'No valid data' message on empty input."""
    x, y = _prep_xy(x_series, y_series)
    if x.size == 0:
        ax.set_title(title, fontsize=16, weight="bold")
        ax.text(0.5, 0.5, "No valid data", ha="center", va="center",
                transform=ax.transAxes, fontsize=16)
        ax.set_xlabel(x_label, fontsize=18, labelpad=10)
        ax.set_ylabel(y_label, fontsize=18, labelpad=10)
        return None, None
    return _scatter_panel(ax, x_series, y_series, title, x_label, y_label)


def _shared_row_limits(panels, pad=0.05):
    """Compute shared axis limits across a list of (x, y) panel results."""
    vals = []
    for x, y in panels:
        if x is not None and y is not None and x.size and y.size:
            vals.extend([x.min(), x.max(), y.min(), y.max()])
    if not vals:
        return None
    lo, hi = float(min(vals)), float(max(vals))
    span   = hi - lo if np.isfinite(hi - lo) else 1.0
    return lo - pad * span, hi + pad * span


def add_esm2(df, index_to_esm2):
    df['esm2'] = pd.array([None] * len(df), dtype=object)
    for ind in df.index:
        esm1bs = []
        for indexes in df.loc[ind, 'index']:
            #print(indexes, type(indexes))
            esm1bs.append(index_to_esm2[str(indexes)])
        df.at[ind, 'esm2'] = esm1bs       
    return df

def add_esm2_max(df, index_to_esm2_max):
    df['esm2_max_pooling'] = pd.array([None] * len(df), dtype=object)
    for ind in df.index:
        esm2s_max = []
        for indexes in df.loc[ind, 'index']:
            #print(indexes, type(indexes))
            esm2s_max.append(index_to_esm2_max[str(indexes)])
        df.at[ind, 'esm2_max_pooling'] = esm2s_max       
    return df

def add_esm2_absmax(df, index_to_esm2_absmax):
    df['esm2_absmax_pooling'] = pd.array([None] * len(df), dtype=object)
    for ind in df.index:
        esm2s_absmax = []
        for indexes in df.loc[ind, 'index']:
            #print(indexes, type(indexes))
            esm2s_absmax.append(index_to_esm2_absmax[str(indexes)])
        df.at[ind, 'esm2_absmax_pooling'] = esm2s_absmax       
    return df

def add_rxnfp(df, index_to_rxnfp):
    df['rxnfp'] = pd.array([None] * len(df), dtype=object)
    for ind in df.index:
        indexes = df.loc[ind, 'index']
        df.at[ind, 'rxnfp'] = index_to_rxnfp[str(indexes[0])]       
    return df

def add_symetric_pairs(df):
    df_sym = df.copy()
    df_sym["Type"] = df_sym["Type"].apply(lambda x: x[::-1])
    df_sym["index"] = df_sym["index"].apply(lambda x: x[::-1])
    df_sym["log10kcat"] = df_sym["log10kcat"].apply(lambda x: x[::-1])
    return df_sym


def prepro(df):
    df['fp1'] = df['esm2'].apply(lambda x: np.asarray(x[0]) -  np.asarray(x[1]))
    df['fp1_max'] = df['esm2_max_pooling'].apply(lambda x: np.asarray(x[0]) -  np.asarray(x[1]))
    df['delta_log_kcat'] = df['log10kcat'].apply(lambda x: x[0] - x[1])
    df['ESM2_1'] = df['esm2'].apply(lambda x: np.asarray(x[0]))
    df['ESM2_max_1'] = df['esm2_max_pooling'].apply(lambda x: np.asarray(x[1]))
    df['fp1_absmax'] = df['esm2_absmax_pooling'].apply(lambda x: np.asarray(x[0]) -  np.asarray(x[1]))
    df['ESM2_absmax_1'] = df['esm2_absmax_pooling'].apply(lambda x: np.asarray(x[0]))
    df['ESM2_absmax_2'] = df['esm2_absmax_pooling'].apply(lambda x: np.asarray(x[1]))
    df = df.dropna(subset=['delta_log_kcat'])
    df = df[~df['delta_log_kcat'].isin([np.inf, -np.inf])]
    df = df.reset_index(drop=True)
    return df

def prepro_v2(df):
    """Preprocessing for the concatenated-embedding model (FCKcat): adds delta_log_kcat and sequences_representation."""
    df = df.copy()
    df['delta_log_kcat'] = df['log10kcat'].apply(lambda x: x[0] - x[1])
    df['sequences_representation'] = df['esm2_absmax_pooling'].apply(
        lambda x: np.concatenate([np.asarray(x[0]), np.asarray(x[1])])
    )
    df['delta_representation'] = df['esm2_absmax_pooling'].apply(
        lambda x: np.asarray(x[0]) - np.asarray(x[1])
    )
    df = df.dropna(subset=['delta_log_kcat'])
    df = df[~df['delta_log_kcat'].isin([np.inf, -np.inf])]
    df = df.reset_index(drop=True)
    return df


def add_esm1v(df, dict_esm1v_max_pooling):
    # create a fresh column
    df['ESM1v_max_pooling'] = pd.array([None] * len(df), dtype=object)

    valid_rows = []  # store rows without wrong_id

    for ind in df.index:
        esm1v_max = []
        wrong = False

        for indexes in df.loc[ind, 'index']:
            value = dict_esm1v_max_pooling.get(str(indexes), None)
            esm1v_max.append(value)
            if isinstance(value, str):  # mark as wrong
                wrong = True

        # only keep row if no wrong_id was found
        if not wrong:
            df.at[ind, 'ESM1v_max_pooling'] = esm1v_max
            valid_rows.append(ind)

    # return df filtered to valid rows only
    return df.loc[valid_rows].copy()


def _shared_row_limits_v2(panels, pad=0.05):
    vals = []
    for x, y in panels:
        if x is not None and y is not None and x.size and y.size:
            vals.extend([x.min(), x.max(), y.min(), y.max()])
    if not vals:
        return None
    lo, hi = float(min(vals)), float(max(vals))
    span = hi - lo if np.isfinite(hi - lo) else 1.0
    return lo - pad*span, hi + pad*span

def _kcat_formatter_from_log10(val, pos):
    """Format log10(kcat) axis values as powers of 10."""
    if not np.isfinite(val):
        return ""

    exponent = int(np.round(val))

    # Only label near-integer ticks to avoid clutter
    if np.isclose(val, exponent, atol=1e-6):
        if exponent == 0:
            return r"$10^{0}$"
        return rf"$10^{{{exponent}}}$"
    else:
        return ""


def _metrics_no_threshold(true_series, pred_series):
    """Compute Pearson r, sign accuracy, and MCC between two series."""
    x, y = _prep_xy(true_series, pred_series)
    n = x.size
    if n < 2:
        return {"pearson": np.nan, "acc": np.nan, "mcc": np.nan, "n": int(n)}

    try:
        r, _ = pearsonr(x, y)
    except Exception:
        r = np.nan

    tp = (x > 0) & (y > 0)
    tn = (x < 0) & (y < 0)
    acc = (tp | tn).sum() / n

    y_true_bin = (x > 0).astype(int)
    y_pred_bin = (y > 0).astype(int)

    try:
        mcc = matthews_corrcoef(y_true_bin, y_pred_bin)
    except Exception:
        mcc = np.nan

    return {
        "pearson": float(r),
        "acc": float(acc),
        "mcc": float(mcc),
        "n": int(n)
    }


def _bars_panel(ax, df, true_col, pred_cols, labels, title=""):
    """Grouped bar chart of Pearson r / Accuracy / MCC for two prediction columns."""
    if title:
        ax.set_title(title, fontsize=20, weight="bold", pad=12)

    res = []
    for pred, lab in zip(pred_cols, labels):
        m = _metrics_no_threshold(df[true_col], df[pred])
        res.append((lab, m))

    metrics_names = ['Pearson r', 'Accuracy', 'MCC']
    vals = {
        res[0][0]: [res[0][1]['pearson'], res[0][1]['acc'], res[0][1]['mcc']],
        res[1][0]: [res[1][1]['pearson'], res[1][1]['acc'], res[1][1]['mcc']],
    }

    x = np.arange(len(metrics_names))
    width = 0.36

    b1 = ax.bar(x - width/2, vals[res[0][0]], width, label=res[0][0])
    b2 = ax.bar(x + width/2, vals[res[1][0]], width, label=res[1][0])

    ax.set_xticks(x)
    ax.set_xticklabels(metrics_names, fontsize=18)
    ax.set_ylabel('Score', fontsize=20, labelpad=10)
    ax.legend(frameon=False, fontsize=18, loc='upper left')

    all_vals = np.array(vals[res[0][0]] + vals[res[1][0]], dtype=float)
    finite_vals = all_vals[np.isfinite(all_vals)]
    ymin = min(-0.1, np.min(finite_vals) - 0.08) if finite_vals.size else -0.1
    ymax = max(1.05, np.max(finite_vals) + 0.20) if finite_vals.size else 1.05
    ax.set_ylim(ymin, ymax)

    for bars in (b1, b2):
        for b in bars:
            h = b.get_height()
            if np.isfinite(h):
                ax.annotate(
                    f'{h:.3f}',
                    (b.get_x() + b.get_width()/2, h),
                    textcoords="offset points",
                    xytext=(0, 6),
                    ha='center',
                    va='bottom',
                    fontsize=15,
                    fontweight='bold'
                )

    ax.tick_params(axis='both', which='major', labelsize=18, width=1.6, length=7)
    ax.tick_params(axis='both', which='minor', width=1.2, length=4)

    for spine in ax.spines.values():
        spine.set_linewidth(1.6)