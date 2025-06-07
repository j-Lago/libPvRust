


const T_REF: f64 = 298.15; // [K] rated temperature
const S_REF: f64 = 1000.0; // [W/m^2] rated irradiance
const Q_K: f64 = 1.60217663e-19 / 1.38064852e-23;  // [K/eV] Boltzmann constante reciprocal
const C_TO_K: f64 = 273.15;

use std::fmt;
use tracing::{warn, error};

pub static mut SOLVER_CALLS: usize = 0;


#[derive(Clone, PartialEq)]
pub struct PvCellSolver {
    pub max_iter: usize,    // max number of iterations
    pub tol_i: f64, // [A] current tolerance
    pub tol_v: f64, // [V] voltage tolerance
}

impl Default for PvCellSolver {
    fn default() -> Self {
        return PvCellSolver { max_iter: 100, tol_i: 0.001, tol_v: 0.01 };
    }
}

#[derive(Clone, Debug)]
pub struct PvCellState {
    pub gsh: f64,
    pub ra: f64,
    pub i0: f64,
    pub il: f64,
}

pub struct BasicParams {
    pub a_ref: f64,
    pub i_o_ref: f64,
    pub i_l_ref: f64,
    pub r_s: f64,
    pub r_sh_ref: f64,
    pub alpha_sc: f64,
    pub v_oc_ref: f64
}

#[derive(Clone, PartialEq)]
pub struct PvCell {
    pub a_ref: f64,
    pub i_o_ref: f64,
    pub i_l_ref: f64,
    pub r_s: f64,
    pub r_sh_ref: f64,
    pub alpha_sc: f64,
    pub v_oc_ref: f64,
    pub v_bypass: f64, // [V] tensão de bypass
    pub r_bypass: f64, // [Ohm] resistance of bypass diode
    pub eg_ref: f64,   // [eV]  band energy Si: 1.121, CdTe: 1.475
    pub degdt: f64,    // Si: -0.0002677, CdTe: -0.0003 //
    pub shading: f64,
    pub np: u32, // [-]   number of modules in parallel
    pub ns: u32, // [-]   number of modules in series
    pub solver: PvCellSolver,
}

impl PvCell {
    pub(crate) fn default() -> PvCell {
        PvCell {
            a_ref: f64::NAN,
            i_o_ref: f64::NAN,
            i_l_ref: f64::NAN,
            r_s: f64::NAN,
            r_sh_ref: f64::NAN,
            alpha_sc: f64::NAN,
            v_oc_ref: f64::NAN,
            v_bypass: -0.65 * 3.0,
            r_bypass: 0.1,
            eg_ref: 1.121,
            degdt: -0.0002677,
            shading: 0.0,
            np: 1,
            ns: 1,
            solver: PvCellSolver::default()
        }
    }
}

#[allow(dead_code)]
impl PvCell {
    pub fn new(params: &BasicParams) -> Self {
        PvCell {
            a_ref: params.a_ref, i_o_ref: params.i_o_ref, i_l_ref: params.i_l_ref, r_s: params.r_s,
            r_sh_ref: params.r_sh_ref, alpha_sc: params.alpha_sc, v_oc_ref: params.v_oc_ref,
            ..PvCell::default()
        }
    }

    /// builders
    pub fn with_ns(mut self, ns: u32) -> Self{ self.ns = ns; return self; }
    pub fn with_np(mut self, np: u32) -> Self{ self.np = np; return self; }
    pub fn with_shading(mut self, shading: f64) -> Self{ self.shading = shading; return self; }
    pub fn with_solver(mut self, settings: PvCellSolver) -> Self { self.solver = settings; return self; }

    pub fn compute_state(&self, irrad_ef: f64, cell_temp: f64) -> PvCellState {
        let irrad: f64 = irrad_ef * (1.0 - self.shading);
        let tj: f64 = cell_temp + C_TO_K;
        let eg: f64 = self.eg_ref * (1. + self.degdt  * (tj - T_REF));
        let gsh: f64 = irrad / (self.r_sh_ref * S_REF);
        let ra: f64 = T_REF / (self.a_ref * tj);
        let i0: f64 = self.i_o_ref  * (tj / T_REF).powi(3) * (Q_K * (self.eg_ref / T_REF - eg / tj)).exp();
        let il: f64 = (self.i_l_ref + self.alpha_sc * (tj - T_REF)) * irrad / S_REF;
        PvCellState {gsh, ra, i0, il}
    }

    pub fn solve_i(&self, state: &PvCellState, v_pnl: f64) -> f64 {
        unsafe{ SOLVER_CALLS += 1}
        let mut i: f64 = 0.0;
        let v: f64 = v_pnl / (self.ns as f64);

        let mut success: bool = false;
        for _ in 0..self.solver.max_iter {
            let den: f64 = -1.0 - state.i0 * ((v + i * self.r_s) * state.ra).exp() * self.r_s * state.ra - self.r_s * state.gsh;
            let d: f64 = (state.il - i - state.i0 * (((v + i * self.r_s) * state.ra).exp() - 1.0) - (v + i * self.r_s) * state.gsh) / den;
            i = i - d;
            if d.abs() < self.solver.tol_i {
                success = true;
                break;
            }
        }
        if v < self.v_bypass {
            i += (self.v_bypass - v) / self.r_bypass;
        }
        i *= self.np as f64;

        if !success {
            if i.is_normal() {
                warn!("({:p}) PvCell::solve_i(v_pnl={:e}) nao convergiu (tol={:e}, max_iter={}) -> (i={})",
                    &self, v_pnl, self.solver.tol_i, self.solver.max_iter, i);
            }else{
                error!("({:p}) PvCell::solve_i(v_pnl={:e}) nao convergiu (tol={:e}, max_iter={}) -> (i={})",
                    &self, v_pnl, self.solver.tol_i, self.solver.max_iter, i);
            }
        }
        return i;
    }

    pub fn v_from_i(&self, state: &PvCellState, i_pnl: f64) -> f64 {
        unsafe{ SOLVER_CALLS += 1}
        let mut v: f64 = self.v_oc_ref;
        let i: f64 = i_pnl / (self.np as f64);

        let mut success = false;
        for _ in 0..self.solver.max_iter {
            let den: f64 = -state.i0 * ((v + i * self.r_s) * state.ra).exp() * state.ra - (state.gsh);
            let d: f64 = (state.il - i - state.i0 * (((v + i * self.r_s) * state.ra).exp() - 1.) - (v + i * self.r_s) * state.gsh) / den;
            v -= d;
            if d.abs() < self.solver.tol_v {
                success = true;
                break;
            }
        }
        if v < self.v_bypass {
            let ir: f64 = self.solve_i(state, self.v_bypass) / (self.np as f64);   // corrente em que se inicia região de breakdown para solve_v(i)
            v = self.v_bypass - (i - ir) * self.r_bypass;
        }
        v *=  self.ns as f64;

        if !success {
            if v.is_normal() {
                warn!("({:p}) PvCell::v_from_i(i_pnl={:e}) nao convergiu (tol={:e}, max_iter={}) -> (v={})",
                    &self, i_pnl, self.solver.tol_v, self.solver.max_iter, v);
            }else{
                error!("({:p}) PvCell::v_from_i(i_pnl={:e}) nao convergiu (tol={:e}, max_iter={}) -> (v={})",
                    &self, i_pnl, self.solver.tol_v, self.solver.max_iter, v);
            }
        }
        return v;
    }

    pub fn is_extended_params_equivalent(&self, other: &PvCell) -> bool {
        self.a_ref == other.a_ref && 
        self.i_o_ref == other.i_o_ref && 
        self.i_l_ref == other.i_l_ref && 
        self.r_s == other.r_s && 
        self.r_sh_ref == other.r_sh_ref && 
        self.alpha_sc == other.alpha_sc && 
        self.v_oc_ref == other.v_oc_ref && 
        self.v_bypass == other.v_bypass && 
        self.r_bypass == other.r_bypass && 
        self.eg_ref == other.eg_ref && 
        self.degdt == other.degdt && 
        self.shading == other.shading
    }
    
    pub fn is_series_equivalent(&self, other: &PvCell) -> bool {
            self.is_extended_params_equivalent(other) && self.np == other.np
    }

    pub fn is_parallel_equivalent(&self, other: &PvCell) -> bool {
        self.is_extended_params_equivalent(other) && self.ns == other.ns
    }
}


impl fmt::Debug for PvCell {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{{{},{}}}", self.ns, self.np)
    }
}