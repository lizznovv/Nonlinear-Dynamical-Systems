// Dynamic_Systems_Bifurcation.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <string>
#include <vector>
#include <cmath> 
#include <fstream>
#include <omp.h>
#include "System.h"
#include "Methods.h"
#include "Bifurcation.h"
#include "Lyapunov.h"
#include "PhaseVolume.h"


int main()
{
    std::vector<double> coefs = { 0.2, 0.2, 5.7, 0.5 };

    /*
    
    std::cout << "\nMethod: Euler";
    std::cout << "\na parameter processing...";
    auto bif_a_euler = build_bifurcation(0.0, 0.36, 1000, 0, coefs, euler);
    save_csv("euler_a.csv", bif_a_euler);
    std::cout << "\ndone";

    std::cout << "\nb parameter processing...";
    auto bif_b_euler = build_bifurcation(0.1, 2.0, 1000, 1, coefs, euler);
    save_csv("euler_b.csv", bif_b_euler);
    std::cout << "\ndone";

    std::cout << "\nc parameter processing...";
    auto bif_c_euler = build_bifurcation(2.0, 10.0, 1000, 2, coefs, euler);
    save_csv("euler_c.csv", bif_c_euler);
    std::cout << "\ndone\n";


    std::cout << "\nMethod: Euler-Cromer";
    std::cout << "\na parameter processing...";
    auto bif_a_ecromer = build_bifurcation(0.0, 0.36, 1000, 0, coefs, e_cromer);
    save_csv("e_cromer_a.csv", bif_a_ecromer);
    std::cout << "\ndone";

    std::cout << "\nb parameter processing...";
    auto bif_b_ecromer = build_bifurcation(0.1, 2.0, 1000, 1, coefs, e_cromer);
    save_csv("e_cromer_b.csv", bif_b_ecromer);
    std::cout << "\ndone";

    std::cout << "\nc parameter processing...";
    auto bif_c_ecromer = build_bifurcation(2.0, 10.0, 1000, 2, coefs, e_cromer);
    save_csv("e_cromer_c.csv", bif_c_ecromer);
    std::cout << "\ndone\n";


    std::cout << "\nMethod: Midpoint";
    std::cout << "\na parameter processing...";
    auto bif_a_midpoint = build_bifurcation(0.0, 0.36, 1000, 0, coefs, midpoint);
    save_csv("midpoint_a.csv", bif_a_midpoint);
    std::cout << "\ndone";

    std::cout << "\nb parameter processing...";
    auto bif_b_midpoint = build_bifurcation(0.1, 2.0, 1000, 1, coefs, midpoint);
    save_csv("midpoint_b.csv", bif_b_midpoint);
    std::cout << "\ndone";

    std::cout << "\nc parameter processing...";
    auto bif_c_midpoint = build_bifurcation(2.0, 10.0, 1000, 2, coefs, midpoint);
    save_csv("midpoint_c.csv", bif_c_midpoint);
    std::cout << "\ndone\n";


    std::cout << "\nMethod: VSCD";
    std::cout << "\na parameter processing...";
    auto bif_a_vscd = build_bifurcation(0.0, 0.36, 1000, 0, coefs, vscd);
    save_csv("vscd_a.csv", bif_a_vscd);
    std::cout << "\ndone";

    std::cout << "\nb parameter processing...";
    auto bif_b_vscd = build_bifurcation(0.1, 2.0, 1000, 1, coefs, vscd);
    save_csv("vscd_b.csv", bif_b_vscd);
    std::cout << "\ndone";

    std::cout << "\nc parameter processing...";
    auto bif_c_vscd = build_bifurcation(2.0, 10.0, 1000, 2, coefs, vscd);
    save_csv("vscd_c.csv", bif_c_vscd);
    std::cout << "\ndone";

    std::cout << "\ns parameter processing...";
    auto bif_s_vscd = build_bifurcation(0.0, 1.0, 1000, 3, coefs, vscd);
    save_csv("vscd_s.csv", bif_s_vscd);
    std::cout << "\ndone\n";


    std::cout << "\nMethod: RK4";
    std::cout << "\na parameter processing...";
    auto bif_a_rk4 = build_bifurcation(0.0, 0.36, 1000, 0, coefs, rk4);
    save_csv("rk4_a.csv", bif_a_rk4);
    std::cout << "\ndone";

    std::cout << "\nb parameter processing...";
    auto bif_b_rk4 = build_bifurcation(0.1, 2.0, 1000, 1, coefs, rk4);
    save_csv("rk4_b.csv", bif_b_rk4);
    std::cout << "\ndone";

    std::cout << "\nc parameter processing...";
    auto bif_c_rk4 = build_bifurcation(2.0, 10.0, 1000, 2, coefs, rk4);
    save_csv("rk4_c.csv", bif_c_rk4);
    std::cout << "\ndone\n";
    
    */

/*
    std::cout << "\nMethod: Euler";
    std::cout << "\nh exploring...";
    auto bif_h_euler = h_build_bifurcation(1e-3, 1e-1, 300, coefs, euler);
    save_csv("euler_h.csv", bif_h_euler);
    std::cout << "\ndone";

    std::cout << "\nMethod: Euler-Cromer";
    std::cout << "\nh exploring...";
    auto bif_h_ecromer = h_build_bifurcation(1e-4, 5e-2, 300, coefs, e_cromer);
    save_csv("e_cromer_h.csv", bif_h_ecromer);
    std::cout << "\ndone";

    std::cout << "\nMethod: Midpoint";
    std::cout << "\nh exploring...";
    auto bif_h_midpoint = h_build_bifurcation(1e-4, 5e-2, 300, coefs, midpoint);
    save_csv("midpoint_h.csv", bif_h_midpoint);
    std::cout << "\ndone";

    std::cout << "\nMethod: VSCD";
    std::cout << "\nh exploring...";
    auto bif_h_vscd = h_build_bifurcation(1e-4, 5e-2, 300, coefs, vscd);
    save_csv("vscd_h.csv", bif_h_vscd);
    std::cout << "\ndone";

    std::cout << "\nMethod: RK4";
    std::cout << "\nh exploring...";
    auto bif_h_rk4 = h_build_bifurcation(1e-4, 5e-2, 300, coefs, rk4);
    save_csv("rk4_h.csv", bif_h_rk4);
    std::cout << "\ndone";

    */

/*
    std::cout << "\nMethod: Euler";
    auto euler_bif_ac = ac_build_bifurcation(0.0, 0.36, 2.0, 10.0, 100, 0, 2, coefs, euler);
    save_2d_csv("ac_euler.csv", euler_bif_ac, "a", "c");
    std::cout << "\ndone";

    std::cout << "\nMethod: Euler-Cromer";
    auto ecromer_bif_ac = ac_build_bifurcation(0.0, 0.36, 2.0, 10.0, 100, 0, 2, coefs, e_cromer);
    save_2d_csv("ac_ecromer.csv", ecromer_bif_ac, "a", "c");
    std::cout << "\ndone";
    
    std::cout << "\nMethod: Midpoint";
    auto midpoint_bif_ac = ac_build_bifurcation(0.0, 0.36, 2.0, 10.0, 100, 0, 2, coefs, midpoint);
    save_2d_csv("ac_midpoint.csv", midpoint_bif_ac, "a", "c");
    std::cout << "\ndone";
    
    std::cout << "\nMethod: VSCD";
    auto vscd_bif_ac = ac_build_bifurcation(0.0, 0.36, 2.0, 10.0, 100, 0, 2, coefs, vscd);
    save_2d_csv("ac_vscd.csv", vscd_bif_ac, "a", "c");
    std::cout << "\ndone";

    std::cout << "\nMethod: RK4";
    auto rk4_bif_ac = ac_build_bifurcation(0.0, 0.36, 2.0, 10.0, 100, 0, 2, coefs, rk4);
    save_2d_csv("ac_rk4.csv", rk4_bif_ac, "a", "c");
    std::cout << "\ndone";
    */

/*
    std::cout << "\nMethod: Euler (a vs h)";
    auto euler_bif_ah = ah_build_bifurcation(0.0, 0.36, 0.001, 0.1, 100, 0, coefs, euler);
    save_2d_csv("ah_euler.csv", euler_bif_ah, "a", "h");
    std::cout << " done";

    std::cout << "\nMethod: Euler-Cromer (a vs h)";
    auto ecromer_bif_ah = ah_build_bifurcation(0.0, 0.36, 0.001, 0.1, 100, 0, coefs, e_cromer);
    save_2d_csv("ah_ecromer.csv", ecromer_bif_ah, "a", "h");
    std::cout << " done";

    std::cout << "\nMethod: Midpoint (a vs h)";
    auto midpoint_bif_ah = ah_build_bifurcation(0.0, 0.36, 0.001, 0.1, 100, 0, coefs, midpoint);
    save_2d_csv("ah_midpoint.csv", midpoint_bif_ah, "a", "h");
    std::cout << " done";

    std::cout << "\nMethod: VSCD (a vs h)";
    auto vscd_bif_ah = ah_build_bifurcation(0.0, 0.36, 0.001, 0.1, 100, 0, coefs, vscd);
    save_2d_csv("ah_vscd.csv", vscd_bif_ah, "a", "h");
    std::cout << " done";

    std::cout << "\nMethod: RK4 (a vs h)";
    auto rk4_bif_ah = ah_build_bifurcation(0.0, 0.36, 0.001, 0.1, 100, 0, coefs, rk4);
    save_2d_csv("ah_rk4.csv", rk4_bif_ah, "a", "h");
    std::cout << " done\n";
    */

/*
    std::cout << "\nEuler";
    auto lle_h_euler = build_LLE_diagram(0.001, 0.1, 300, 0, 400, 150, coefs, state, euler, true);
    std::cout << "\nSaving Euler...";
    save_lle_csv("lle_euler_h.csv", lle_h_euler);
    std::cout << " done";

    std::cout << "\nEuler-Cromer";
    auto lle_h_ecromer = build_LLE_diagram(0.001, 0.1, 300, 0, 400, 150, coefs, state, e_cromer, true);
    std::cout << "\nSaving EulerCromer...";
    save_lle_csv("lle_ecromer_h.csv", lle_h_ecromer);
    std::cout << " done";

    std::cout << "\nMidpoint";
    auto lle_h_midpoint = build_LLE_diagram(0.001, 0.1, 300, 0, 400, 150, coefs, state, midpoint, true);
    std::cout << "\nSaving Midpoint...";
    save_lle_csv("lle_midpoint_h.csv", lle_h_midpoint);
    std::cout << " done";

    std::cout << "\nVSCD";
    auto lle_h_vscd = build_LLE_diagram(0.001, 0.1, 300, 0, 400, 150, coefs, state, vscd, true);
    std::cout << "\nSaving VSCD...";
    save_lle_csv("lle_vscd_h.csv", lle_h_vscd);
    std::cout << " done";

    std::cout << "\nRK4";
    auto lle_h_rk4 = build_LLE_diagram(0.001, 0.1, 300, 0, 400, 150, coefs, state, rk4, true);
    std::cout << "\nSaving RK4...";
    save_lle_csv("lle_rk4_h.csv", lle_h_rk4);
    std::cout << " done";

    */

/*
    auto ac_lle_euler = ac_build_LLE_diagram(0.0, 0.36, 2.0, 10.0, 100, 0, 2, 400, 150, coefs, state, euler);
    save_lle_map_csv("ac_lle_euler.csv", ac_lle_euler);

    auto ac_lle_ceuler = ac_build_LLE_diagram(0.0, 0.36, 2.0, 10.0, 100, 0, 2, 400, 150, coefs, state, e_cromer);
    save_lle_map_csv("ac_lle_ceuler.csv", ac_lle_ceuler);

    auto ac_lle_midpoint = ac_build_LLE_diagram(0.0, 0.36, 2.0, 10.0, 100, 0, 2, 400, 150, coefs, state, midpoint);
    save_lle_map_csv("ac_lle_midpoint.csv", ac_lle_midpoint);

    auto ac_lle_vscd = ac_build_LLE_diagram(0.0, 0.36, 2.0, 10.0, 100, 0, 2, 400, 150, coefs, state, vscd);
    save_lle_map_csv("ac_lle_vscd.csv", ac_lle_vscd);

    auto ac_lle_rk4 = ac_build_LLE_diagram(0.0, 0.36, 2.0, 10.0, 100, 0, 2, 400, 150, coefs, state, rk4);
    save_lle_map_csv("ac_lle_rk4.csv", ac_lle_rk4);
    */

/*
    auto ah_lle_euler = ah_build_LLE_diagram(0.0, 0.36, 0.001, 0.02, 50, 0, 400, 150, coefs, state, euler);
    save_lle_map_csv("ah_lle_euler.csv", ah_lle_euler);

    auto ah_lle_ecromer = ah_build_LLE_diagram(0.0, 0.36, 0.001, 0.02, 50, 0, 400, 150, coefs, state, e_cromer);
    save_lle_map_csv("ah_lle_ecromer.csv", ah_lle_ecromer);

    auto ah_lle_midpoint = ah_build_LLE_diagram(0.0, 0.36, 0.001, 0.02, 50, 0, 400, 150, coefs, state, midpoint);
    save_lle_map_csv("ah_lle_midpoint.csv", ah_lle_midpoint);

    auto ah_lle_vscd = ah_build_LLE_diagram(0.0, 0.36, 0.001, 0.02, 50, 0, 400, 150, coefs, state, vscd);
    save_lle_map_csv("ah_lle_vscd.csv", ah_lle_vscd);

    auto ah_lle_rk4 = ah_build_LLE_diagram(0.0, 0.36, 0.001, 0.02, 50, 0, 400, 150, coefs, state, rk4);
    save_lle_map_csv("ah_lle_rk4.csv", ah_lle_rk4);
*/

    /*
    std::cout << "\nMethod: Euler";
    std::cout << "\na parameter processing...";
    auto volume_a_euler = build_volume_diagram(0.0, 0.36, 300, 0, 400, 150, coefs, state, euler, false);
    save_volume_csv("volume_euler_a.csv", volume_a_euler);
    std::cout << "\ndone";

    std::cout << "\nb parameter processing...";
    auto volume_b_euler = build_volume_diagram(0.1, 2.0, 300, 1, 400, 150, coefs, state, euler, false);
    save_volume_csv("volume_euler_b.csv", volume_b_euler);
    std::cout << "\ndone";

    std::cout << "\nc parameter processing...";
    auto volume_c_euler = build_volume_diagram(2.0, 10.0, 300, 2, 400, 150, coefs, state, euler, false);
    save_volume_csv("volume_euler_c.csv", volume_c_euler);
    std::cout << "\ndone";

    std::cout << "\nh parameter processing...";
    auto volume_h_euler = build_volume_diagram(0.001, 0.02, 300, 0, 400, 150, coefs, state, euler, true);
    save_volume_csv("volume_euler_h.csv", volume_h_euler);
    std::cout << "\ndone\n";


    std::cout << "\nMethod: Euler-Cromer";
    std::cout << "\na parameter processing...";
    auto volume_a_ecromer = build_volume_diagram(0.0, 0.36, 300, 0, 400, 150, coefs, state, e_cromer, false);
    save_volume_csv("volume_ecromer_a.csv", volume_a_ecromer);
    std::cout << "\ndone";

    std::cout << "\nb parameter processing...";
    auto volume_b_ecromer = build_volume_diagram(0.1, 2.0, 300, 1, 400, 150, coefs, state, e_cromer, false);
    save_volume_csv("volume_ecromer_b.csv", volume_b_ecromer);
    std::cout << "\ndone";

    std::cout << "\nc parameter processing...";
    auto volume_c_ecromer = build_volume_diagram(2.0, 10.0, 300, 2, 400, 150, coefs, state, e_cromer, false);
    save_volume_csv("volume_ecromer_c.csv", volume_c_ecromer);
    std::cout << "\ndone";

    std::cout << "\nh parameter processing...";
    auto volume_h_ecromer = build_volume_diagram(0.001, 0.02, 300, 0, 400, 150, coefs, state, e_cromer, true);
    save_volume_csv("volume_ecromer_h.csv", volume_h_ecromer);
    std::cout << "\ndone\n";


    std::cout << "\nMethod: Midpoint";
    std::cout << "\na parameter processing...";
    auto volume_a_midpoint = build_volume_diagram(0.0, 0.36, 300, 0, 400, 150, coefs, state, midpoint, false);
    save_volume_csv("volume_midpoint_a.csv", volume_a_midpoint);
    std::cout << "\ndone";

    std::cout << "\nb parameter processing...";
    auto volume_b_midpoint = build_volume_diagram(0.1, 2.0, 300, 1, 400, 150, coefs, state, midpoint, false);
    save_volume_csv("volume_midpoint_b.csv", volume_b_midpoint);
    std::cout << "\ndone";

    std::cout << "\nc parameter processing...";
    auto volume_c_midpoint = build_volume_diagram(2.0, 10.0, 300, 2, 400, 150, coefs, state, midpoint, false);
    save_volume_csv("volume_midpoint_c.csv", volume_c_midpoint);
    std::cout << "\ndone";

    std::cout << "\nh parameter processing...";
    auto volume_h_midpoint = build_volume_diagram(0.001, 0.02, 300, 0, 400, 150, coefs, state, midpoint, true);
    save_volume_csv("volume_midpoint_h.csv", volume_h_midpoint);
    std::cout << "\ndone\n";


    std::cout << "\nMethod: VSCD";
    std::cout << "\na parameter processing...";
    auto volume_a_vscd = build_volume_diagram(0.0, 0.36, 300, 0, 400, 150, coefs, state, vscd, false);
    save_volume_csv("volume_vscd_a.csv", volume_a_vscd);
    std::cout << "\ndone";

    std::cout << "\nb parameter processing...";
    auto volume_b_vscd = build_volume_diagram(0.1, 2.0, 300, 1, 400, 150, coefs, state, vscd, false);
    save_volume_csv("volume_vscd_b.csv", volume_b_vscd);
    std::cout << "\ndone";

    std::cout << "\nc parameter processing...";
    auto volume_c_vscd = build_volume_diagram(2.0, 10.0, 300, 2, 400, 150, coefs, state, vscd, false);
    save_volume_csv("volume_vscd_c.csv", volume_c_vscd);
    std::cout << "\ndone";

    std::cout << "\ns parameter processing...";
    auto volume_s_vscd = build_volume_diagram(0.0, 1.0, 300, 3, 400, 150, coefs, state, vscd, false);
    save_volume_csv("volume_vscd_s.csv", volume_s_vscd);
    std::cout << "\ndone";

    std::cout << "\nh parameter processing...";
    auto volume_h_vscd = build_volume_diagram(0.001, 0.02, 300, 0, 400, 150, coefs, state, vscd, true);
    save_volume_csv("volume_vscd_h.csv", volume_h_vscd);
    std::cout << "\ndone\n";


    std::cout << "\nMethod: RK4";
    std::cout << "\na parameter processing...";
    auto volume_a_rk4 = build_volume_diagram(0.0, 0.36, 300, 0, 400, 150, coefs, state, rk4, false);
    save_volume_csv("volume_rk4_a.csv", volume_a_rk4);
    std::cout << "\ndone";

    std::cout << "\nb parameter processing...";
    auto volume_b_rk4 = build_volume_diagram(0.1, 2.0, 300, 1, 400, 150, coefs, state, rk4, false);
    save_volume_csv("volume_rk4_b.csv", volume_b_rk4);
    std::cout << "\ndone";

    std::cout << "\nc parameter processing...";
    auto volume_c_rk4 = build_volume_diagram(2.0, 10.0, 300, 2, 400, 150, coefs, state, rk4, false);
    save_volume_csv("volume_rk4_c.csv", volume_c_rk4);
    std::cout << "\ndone";

    std::cout << "\nh parameter processing...";
    auto volume_h_rk4 = build_volume_diagram(0.001, 0.02, 300, 0, 400, 150, coefs, state, rk4, true);
    save_volume_csv("volume_rk4_h.csv", volume_h_rk4);
    std::cout << "\ndone\n";
  */

    std::cout << "\nMethod: Euler";

    std::cout << "\nh = 0.1...";
    auto cont_euler_01 = build_contunue_diagram(0.1, 1000, state, 1000000, 10000, euler, coefs);
    save_csv("continue_euler_h01.csv", cont_euler_01);
    std::cout << "\ndone";

    std::cout << "\nh = 0.05...";
    auto cont_euler_005 = build_contunue_diagram(0.05, 1000, state, 1000000, 10000, euler, coefs);
    save_csv("continue_euler_h005.csv", cont_euler_005);
    std::cout << "\ndone";

    std::cout << "\nh = 0.01...";
    auto cont_euler_001 = build_contunue_diagram(0.01, 1000, state, 1000000, 10000, euler, coefs);
    save_csv("continue_euler_h001.csv", cont_euler_001);
    std::cout << "\ndone\n";


    std::cout << "\nMethod: Euler-Cromer";

    std::cout << "\nh = 0.1...";
    auto cont_ecromer_01 = build_contunue_diagram(0.1, 1000, state, 1000000, 10000, e_cromer, coefs);
    save_csv("continue_ecromer_h01.csv", cont_ecromer_01);
    std::cout << "\ndone";

    std::cout << "\nh = 0.05...";
    auto cont_ecromer_005 = build_contunue_diagram(0.05, 1000, state, 1000000, 10000, e_cromer, coefs);
    save_csv("continue_ecromer_h005.csv", cont_ecromer_005);
    std::cout << "\ndone";

    std::cout << "\nh = 0.01...";
    auto cont_ecromer_001 = build_contunue_diagram(0.01, 1000, state, 1000000, 10000, e_cromer, coefs);
    save_csv("continue_ecromer_h001.csv", cont_ecromer_001);
    std::cout << "\ndone\n";


    std::cout << "\nMethod: Midpoint";

    std::cout << "\nh = 0.1...";
    auto cont_midpoint_01 = build_contunue_diagram(0.1, 1000, state, 1000000, 10000, midpoint, coefs);
    save_csv("continue_midpoint_h01.csv", cont_midpoint_01);
    std::cout << "\ndone";

    std::cout << "\nh = 0.05...";
    auto cont_midpoint_005 = build_contunue_diagram(0.05, 1000, state, 1000000, 10000, midpoint, coefs);
    save_csv("continue_midpoint_h005.csv", cont_midpoint_005);
    std::cout << "\ndone";

    std::cout << "\nh = 0.01...";
    auto cont_midpoint_001 = build_contunue_diagram(0.01, 1000, state, 1000000, 10000, midpoint, coefs);
    save_csv("continue_midpoint_h001.csv", cont_midpoint_001);
    std::cout << "\ndone\n";


    std::cout << "\nMethod: VSCD";

    std::cout << "\nh = 0.1...";
    auto cont_vscd_01 = build_contunue_diagram(0.1, 1000, state, 1000000, 10000, vscd, coefs);
    save_csv("continue_vscd_h01.csv", cont_vscd_01);
    std::cout << "\ndone";

    std::cout << "\nh = 0.05...";
    auto cont_vscd_005 = build_contunue_diagram(0.05, 1000, state, 1000000, 10000, vscd, coefs);
    save_csv("continue_vscd_h005.csv", cont_vscd_005);
    std::cout << "\ndone";

    std::cout << "\nh = 0.01...";
    auto cont_vscd_001 = build_contunue_diagram(0.01, 1000, state, 1000000, 10000, vscd, coefs);
    save_csv("continue_vscd_h001.csv", cont_vscd_001);
    std::cout << "\ndone\n";


    std::cout << "\nMethod: RK4";

    std::cout << "\nh = 0.1...";
    auto cont_rk4_01 = build_contunue_diagram(0.1, 1000, state, 1000000, 10000, rk4, coefs);
    save_csv("continue_rk4_h01.csv", cont_rk4_01);
    std::cout << "\ndone";

    std::cout << "\nh = 0.05...";
    auto cont_rk4_005 = build_contunue_diagram(0.05, 1000, state, 1000000, 10000, rk4, coefs);
    save_csv("continue_rk4_h005.csv", cont_rk4_005);
    std::cout << "\ndone";

    std::cout << "\nh = 0.01...";
    auto cont_rk4_001 = build_contunue_diagram(0.01, 1000, state, 1000000, 10000, rk4, coefs);
    save_csv("continue_rk4_h001.csv", cont_rk4_001);
    std::cout << "\ndone\n";


    return 0;
}

