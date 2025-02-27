import numpy as np
from numpy import float64
from numpy.typing import NDArray


def sfdata(ts:NDArray[float64],
           vals:NDArray[float64],
           errs:NDArray[float64],
           z:float64=0):
    mask = ~(np.isnan(ts) | np.isnan(vals) | np.isnan(errs))
    ts = ts[mask]
    vals = vals[mask]
    errs = errs[mask]


    if len(ts) < 2:
        return np.array([]), np.array([]), np.array([])
    
    taus = ts[1:] - ts[0]

    dvals = np.abs(vals[1:] - vals[0])
    dval_errs = np.sqrt(
        np.square(errs[1:])
            +
        np.square(errs[0])
    )

    for i in range(1, len(ts)-1):
        taus = np.concatenate([taus, ts[i+1:] - ts[i]])
        dvals = np.concatenate([dvals, np.abs(vals[i+1:] - vals[i])])
        dval_errs = np.concatenate([
            dval_errs,
            np.sqrt(
                np.square(errs[i+1:])
                    +
                np.square(errs[i])
            )
        ])
    
    # mask = ~(np.isnan(dvals) | np.isnan(dval_errs))
    # return taus[mask]/(1+z), dvals[mask], dval_errs[mask]
    return taus/(1+z), dvals, dval_errs


def esfdata(time_series_list:list[dict]):
    taus      = np.array([], dtype=np.float64)
    dvals     = np.array([], dtype=np.float64)
    dval_errs = np.array([], dtype=np.float64)

    for series in time_series_list:
        redshift = series.get("z", 0)
        ts = series.get("ts", np.array([]))
        vals = series.get("vals", np.array([]))
        errs = series.get("errs", np.array([]))

        _taus, _dvals, _dval_errs = sfdata(ts, vals, errs, redshift)
        taus      = np.concatenate([taus,      _taus])
        dvals     = np.concatenate([dvals,     _dvals])
        dval_errs = np.concatenate([dval_errs, _dval_errs])
    
    return taus, dvals, dval_errs


def calc_sf(taus:NDArray[float64], dvals:NDArray[float64], dval_errs:NDArray[float64],
            step:float, factor:float=np.pi/2, start:float=0):

    num_taubin = int(np.ceil(np.max(taus)/step))

    X     = np.empty(num_taubin, dtype=np.float64)
    V     = np.empty(num_taubin, dtype=np.float64)
    V_err = np.empty(num_taubin, dtype=np.float64)

    for i in range(0, num_taubin):
        mask = (taus > step*i+start) & (taus <= step*(i+1)+start)
        _dvals = dvals[mask]
        _dval_errs = dval_errs[mask]
        _taus = taus[mask]

        X[i] = np.median(_taus)

        if len(_dvals) > 1:
            mean_dval = np.mean(_dvals)
            
            expr1 = factor * mean_dval * mean_dval
            expr2 = np.mean(_dval_errs**2)

            if expr1 >= expr2:
                V[i] = np.sqrt(expr1-expr2)
                sig_sq_1 = np.var(_dvals, ddof=1)/len(_dvals)
                sig_sq_2 = np.var(dval_errs**2, ddof=1)/len(dval_errs)
                V_err[i] = np.sqrt((np.pi**2) * (mean_dval**2) * sig_sq_1 + sig_sq_2)/(2*V[i])
            else:
                V[i] = np.nan
                V_err[i] = np.nan
        else:
            V[i] = np.nan
            V_err[i] = np.nan

    return X, V, V_err

def wise_sf(time_series:dict, factor=np.pi/2):
    ts = time_series.get("ts")
    vals = time_series.get("vals")
    errs = time_series.get("errs")
    z = time_series.get("z", 0)

    taus, dvals, dval_errs = sfdata(ts, vals, errs, z)

    X = np.unique(taus)
    X = np.sort(X)

    V     = np.empty(len(X), dtype=np.float64)
    for i in range(0, len(X)):
        mask = taus == X[i]
        _dvals = dvals[mask]
        _dval_errs = dval_errs[mask]

        mean_dval = np.mean(_dvals)
        expr1 = factor * (mean_dval*mean_dval)
        expr2 = np.mean(_dval_errs**2)

        if expr1 >= expr2:
            V[i] = np.sqrt(expr1-expr2)
        else:
            V[i] = np.nan

    return X, V

def sf(time_series:dict, step:float, factor:float=np.pi/2, start:float=0):
    ts = time_series.get("ts")
    vals = time_series.get("vals")
    errs = time_series.get("errs")
    z = time_series.get("z", 0)

    taus, dvals, dval_errs = sfdata(ts, vals, errs, z)

    X, V, V_err = calc_sf(taus=taus, dvals=dvals, dval_errs=dval_errs,
                                    step=step, factor=factor, start=start)

    mask = (~np.isnan(V)) & (~np.isnan(V_err))

    X = X[mask]
    V = V[mask]
    V_err = V_err[mask]

    return X, V, V_err


def esf(time_series_list:list[dict], step:float, factor:float=np.pi/2, start:float=0):
    taus, dvals, dval_errs = esfdata(time_series_list)

    X, V, V_err = calc_sf(taus=taus, dvals=dvals, dval_errs=dval_errs,
                                    step=step, factor=factor, start=start)

    mask = (~np.isnan(V)) & (~np.isnan(V_err))

    X = X[mask]
    V = V[mask]
    V_err = V_err[mask]

    return X, V, V_err
    