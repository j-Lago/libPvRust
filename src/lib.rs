mod pvcell;
mod series;

#[cfg(test)]
mod tests {
    use super::*;

    use crate::pvcell::*;
    use crate::series::*;
    use std::time::Instant;
    use tracing_subscriber::fmt::writer::BoxMakeWriter;
    use tracing::{info};
    use std::io::{self, BufRead};
    use std::path::Path;
    use tracing::subscriber::with_default;


    #[test]
    fn test_original(){
        let res = run(1, "log_original.txt", false);
        assert!((res - 49475.47731878234).abs() < 1.0);
    }

    #[test]
    fn test_reduced(){
        let res = run(1, "log_reduced.txt", true);
        assert!((res - 49475.47731878234).abs() < 1.0);
    }

    fn run(repeat_n: usize, filename: &str, reduce: bool) -> f64{

        let mut sum_p: f64 = 0.0;

        let params0 = Params::Basic{
            a_ref: 1.81,
            i_o_ref: 8.5e-11,
            i_l_ref: 7.4,
            r_s: 0.6,
            r_sh_ref: 600.0,
            alpha_sc: 3.8e-3,
            v_oc_ref:48.6,
        };
        let params1 = Params::Extended{
            a_ref: 1.94,
            i_o_ref: 2.5e-11,
            i_l_ref: 9.3,
            r_s: 0.4,
            r_sh_ref: 600.0,
            alpha_sc: 3.8e-3,
            v_oc_ref: 47.4,
            v_bypass: -0.65 * 3.0,
            r_bypass: 0.1,
            eg_ref: 1.121,
            degdt: -0.0002677,
        };

        let solver_settings = PvCellSolver {
            max_iter: 100,
            tol_i: 0.001,
            tol_v: 0.01,
        };

        let pnl0: PvCell = PvCell::new(&params0).with_ns(3).with_np(1).with_solver(solver_settings);
        let pnl1: PvCell = PvCell::new(&params1).with_np(2);
        let pnl2: PvCell = PvCell::new(&Params::Basic{
            a_ref: 1.94,
            i_o_ref: 2.5e-11,
            i_l_ref: 9.3,
            r_s: 0.4,
            r_sh_ref: 600.0,
            alpha_sc: 3.8e-3,
            v_oc_ref: 47.4,
        });

        let mut string: Series = Series::new(vec![
            pnl0.clone(),
            pnl1.clone(),
            pnl2.clone(),
            pnl2.clone(),
            pnl0.clone(),
            pnl2.clone(),
            pnl1.clone(),
            pnl0.clone(),
            pnl0.clone(),
            pnl0.clone(),])
            .with_solver(SeriesSolver {max_iter: 1000, tol_v: 1e-1, min_g: 0.0});
        println!("\noriginal: {:?}", &string);

        let conditions = [(200.0, 60.0), (800.0, 30.0), (600.0, 25.0), (999.0, 45.0)];
        let voltages = [40.0*6.0, 80.0*6.0, 90.0*6.0, 110.0*6.0, 125.0*6.0];

        let log_path_str = format!("./log/{}", filename);
        let log_path = Path::new(log_path_str.as_str());

        let file = std::fs::File::create(log_path).unwrap();
        let subscriber = tracing_subscriber::fmt()
            .with_writer(BoxMakeWriter::new(file))
            .finish();
        with_default(subscriber, || {
            // tracing::subscriber::set_global_default(subscriber).expect("Erro ao configurar logger");
        info!("Log iniciado: {}", log_path.display());

        let start = Instant::now();
        if reduce{
            info!("({:p}) Series::reduce()", &string);
            let (reduced, map_o_to_r, map_r_to_o) = string.reduce();
            string = reduced;
            println!(" reduced: {:?}", &string);
            println!("map_o_to_r: {:?}", &map_o_to_r);
            println!("map_r_to_o: {:?}", &map_r_to_o);
        }



        for _ in 0..repeat_n {
            for cond in conditions {
                let (irrad, temp) = cond;
                let states: Vec<CellState> = string.states_uniform_conditions(irrad, temp);
                for v in voltages {
                    let i = string.i_from_v(&states, v);
                    let v2 = string.v_from_i(&states, i);
                    sum_p += v2 * i;
                    // println!("({irrad}, {temp}): iters: {iter}   {v:.3} -> {i:.3} -> {v2:.3},     pot: {:.3}    dv: {:.3}", i * v2, (v - v2).abs())
                }
            }
        }
        let elapsed_time = start.elapsed().as_secs_f64();

        unsafe{
            let string_count = series::SOLVER_CALLS;
            info!("string::SOLVER_CALLS: {}", string_count);

            let panel_count = pvcell::SOLVER_CALLS;
            info!("panel::SOLVER_CALLS: {}", panel_count);
        }
        info!("elapsed time: {} s", &elapsed_time);
        info!("sum_p: {sum_p}");
        info!("Log encerrado");
        let file = std::fs::File::open(log_path).unwrap();
        let reader = io::BufReader::new(file);

        for line in reader.lines() {
            if let Ok(linha) = line {
                println!("{}", linha);
            }
        }
        });
        
        return sum_p;
    }
}