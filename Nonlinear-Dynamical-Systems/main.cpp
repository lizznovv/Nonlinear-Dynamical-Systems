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


int main()
{
    std::vector<double> coefs = { 0.2, 0.2, 5.7, 0.5 };
    /*std::vector<std::string> names = {
        "euler_a.csv", "euler_b.csv", "euler_c.csv",
        "e_cromer_a.csv", "e_cromer_b.csv", "e_cromer_c.csv",
        "midpoint_a.csv", "midpoint_b.csv", "midpoint_c.csv",
        "vscd_a.csv", "vscd_b.csv", "vscd_c.csv", "vscd_s.csv",
        "rk4_a.csv", "rk4_b.csv", "rk4_c.csv",
    };
    */

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

    std::cout << "\nEuler";
    auto lle_a_euler = build_LLE_diagram(0.0, 0.36, 300, 0, 40000, 400, 150, coefs, state, euler);
    std::cout << "\nSaving Euler...";
    save_lle_csv("lle_euler_a.csv", lle_a_euler);
    std::cout << " done";

    std::cout << "\nEuler-Cromer";
    auto lle_a_ecromer = build_LLE_diagram(0.0, 0.36, 300, 0, 40000, 400, 150, coefs, state, e_cromer);
    std::cout << "\nSaving EulerCromer...";
    save_lle_csv("lle_ecromer_a.csv", lle_a_ecromer);
    std::cout << " done";

    std::cout << "\nMidpoint";
    auto lle_a_midpoint = build_LLE_diagram(0.0, 0.36, 300, 0, 40000, 400, 150, coefs, state, midpoint);
    std::cout << "\nSaving Midpoint...";
    save_lle_csv("lle_midpoint_a.csv", lle_a_midpoint);
    std::cout << " done";

    std::cout << "\nVSCD";
    auto lle_a_vscd = build_LLE_diagram(0.0, 0.36, 300, 0, 40000, 400, 150, coefs, state, vscd);
    std::cout << "\nSaving VSCD...";
    save_lle_csv("lle_vscd_a.csv", lle_a_vscd);
    std::cout << " done";

    std::cout << "\nRK4";
    auto lle_a_rk4 = build_LLE_diagram(0.0, 0.36, 300, 0, 40000, 400, 150, coefs, state, rk4);
    std::cout << "\nSaving RK4...";
    save_lle_csv("lle_rk4_a.csv", lle_a_rk4);
    std::cout << " done";

    return 0;
}

