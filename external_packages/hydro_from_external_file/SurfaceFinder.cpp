// Copyright @ Chun Shen 2014
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>
#include <fstream>

#include "./SurfaceFinder.h"
#include <Cornelius/Cornelius.h>
#include <array>
#include <memory>
using JetscapeCornelius::Cornelius;


using namespace std;

SurfaceFinder::SurfaceFinder(void* hydroinfo_ptr_in,
                             ParameterReader* paraRdr_in) {
    paraRdr = paraRdr_in;
    hydro_type = paraRdr->getVal("hydro_type");
    if (hydro_type == 0) {
#ifdef USE_HDF5
        hydroinfo_ptr = (HydroinfoH5*) hydroinfo_ptr_in;
#endif
    } else {
        hydroinfo_MUSIC_ptr = (Hydroinfo_MUSIC*) hydroinfo_ptr_in;
    }
    T_cut = paraRdr->getVal("T_cut");
}

SurfaceFinder::SurfaceFinder(void* hydroinfo_ptr_in,
                             ParameterReader* paraRdr_in, double T_cut_in) {
    paraRdr = paraRdr_in;
    hydro_type = paraRdr->getVal("hydro_type");
    if (hydro_type == 0) {
#ifdef USE_HDF5
        hydroinfo_ptr = (HydroinfoH5*) hydroinfo_ptr_in;
#endif
    } else {
        hydroinfo_MUSIC_ptr = (Hydroinfo_MUSIC*) hydroinfo_ptr_in;
    }
    T_cut = T_cut_in;
}

SurfaceFinder::~SurfaceFinder() {}

bool SurfaceFinder::check_intersect(double T_cut, double tau, double x,
                                    double y, double dt, double dx, double dy,
                                    std::array<std::array<std::array<double, 2>, 2>, 2> &cube) {
    hydrofluidCell *fluidCellptr = new hydrofluidCell();
    bool intersect = true;

    double tau_low = tau - dt/2.;
    double tau_high = tau + dt/2.;
    double x_left = x - dx/2.;
    double x_right = x + dx/2.;
    double y_left = y - dy/2.;
    double y_right = y + dy/2.;

    if (hydro_type == 0) {
#ifdef USE_HDF5
        hydroinfo_ptr->getHydroinfo(tau_low, x_left, y_left, fluidCellptr);
#endif
    } else {
        hydroinfo_MUSIC_ptr->getHydroValues(x_left, y_left, 0.0, tau_low,
                                            fluidCellptr);
    }
    cube[0][0][0] = fluidCellptr->temperature;
    if (hydro_type == 0) {
#ifdef USE_HDF5
        hydroinfo_ptr->getHydroinfo(tau_low, x_left, y_right, fluidCellptr);
#endif
    } else {
        hydroinfo_MUSIC_ptr->getHydroValues(x_left, y_right, 0.0, tau_low,
                                            fluidCellptr);
    }
    cube[0][0][1] = fluidCellptr->temperature;
    if (hydro_type == 0) {
#ifdef USE_HDF5
        hydroinfo_ptr->getHydroinfo(tau_low, x_right, y_left, fluidCellptr);
#endif
    } else {
        hydroinfo_MUSIC_ptr->getHydroValues(x_right, y_left, 0.0, tau_low,
                                            fluidCellptr);
    }
    cube[0][1][0] = fluidCellptr->temperature;
    if (hydro_type == 0) {
#ifdef USE_HDF5
        hydroinfo_ptr->getHydroinfo(tau_low, x_right, y_right, fluidCellptr);
#endif
    } else {
        hydroinfo_MUSIC_ptr->getHydroValues(x_right, y_right, 0.0, tau_low,
                                            fluidCellptr);
    }
    cube[0][1][1] = fluidCellptr->temperature;
    if (hydro_type == 0) {
#ifdef USE_HDF5
        hydroinfo_ptr->getHydroinfo(tau_high, x_left, y_left, fluidCellptr);
#endif
    } else {
        hydroinfo_MUSIC_ptr->getHydroValues(x_left, y_left, 0.0, tau_high,
                                            fluidCellptr);
    }
    cube[1][0][0] = fluidCellptr->temperature;
    if (hydro_type == 0) {
#ifdef USE_HDF5
        hydroinfo_ptr->getHydroinfo(tau_high, x_left, y_right, fluidCellptr);
#endif
    } else {
        hydroinfo_MUSIC_ptr->getHydroValues(x_left, y_right, 0.0, tau_high,
                                            fluidCellptr);
    }
    cube[1][0][1] = fluidCellptr->temperature;
    if (hydro_type == 0) {
#ifdef USE_HDF5
        hydroinfo_ptr->getHydroinfo(tau_high, x_right, y_left, fluidCellptr);
#endif
    } else {
        hydroinfo_MUSIC_ptr->getHydroValues(x_right, y_left, 0.0, tau_high,
                                            fluidCellptr);
    }
    cube[1][1][0] = fluidCellptr->temperature;
    if (hydro_type == 0) {
#ifdef USE_HDF5
        hydroinfo_ptr->getHydroinfo(tau_high, x_right, y_right, fluidCellptr);
#endif
    } else {
        hydroinfo_MUSIC_ptr->getHydroValues(x_right, y_right, 0.0, tau_high,
                                            fluidCellptr);
    }
    cube[1][1][1] = fluidCellptr->temperature;

    if ((T_cut - cube[0][0][0])*(cube[1][1][1] - T_cut) < 0.0)
        if ((T_cut - cube[0][1][0])*(cube[1][0][1] - T_cut) < 0.0)
            if ((T_cut - cube[0][1][1])*(cube[1][0][0] - T_cut) < 0.0)
                if ((T_cut - cube[0][0][1])*(cube[1][1][0] - T_cut) < 0.0)
                    intersect = false;

    delete fluidCellptr;
    return(intersect);
}

int SurfaceFinder::Find_full_hypersurface() {
    ofstream output;
    output.open("hyper_surface_2+1d.dat");

    double grid_tau0 = 0.0;
    double grid_tauf = 0.0;
    double grid_x0 = 0.0;
    double grid_y0 = 0.0;
    if (hydro_type == 1) {
#ifdef USE_HDF5
        grid_tau0 = hydroinfo_ptr->getHydrogridTau0();
        grid_tauf = hydroinfo_ptr->getHydrogridTaumax();
        grid_x0 = hydroinfo_ptr->getHydrogridX0();
        grid_y0 = hydroinfo_ptr->getHydrogridY0();
#endif
    } else {
        grid_tau0 = hydroinfo_MUSIC_ptr->get_hydro_tau0();
        grid_tauf = hydroinfo_MUSIC_ptr->get_hydro_tau_max();
        grid_x0 = (- hydroinfo_MUSIC_ptr->get_hydro_x_max()
                   + hydroinfo_MUSIC_ptr->get_hydro_dx());
        grid_y0 = grid_x0;
    }

    double grid_dt = paraRdr->getVal("grid_dt");
    double grid_dx = paraRdr->getVal("grid_dx");
    double grid_dy = paraRdr->getVal("grid_dy");

    int dim = 3;
    std::array<double, 4> lattice_spacing_arr{{grid_dt, grid_dx, grid_dy, 1.0}};

    std::unique_ptr<Cornelius> cornelius_ptr = std::make_unique<Cornelius>();
    cornelius_ptr->init_cornelius(dim, T_cut, lattice_spacing_arr);
  
    int ntime = static_cast<int>((grid_tauf - grid_tau0)/grid_dt);
    int nx = static_cast<int>(fabs(2.*grid_x0)/grid_dx);
    int ny = static_cast<int>(fabs(2.*grid_y0)/grid_dy);

    std::array<std::array<std::array<double, 2>, 2>, 2> cube{};
    
    hydrofluidCell *fluidCellptr = new hydrofluidCell();
  
    for (int itime = 0; itime < ntime; itime++) {
        // loop over time evolution
        double tau_local = grid_tau0 + (itime + 0.5)*grid_dt;
        for (int i = 0; i < nx; i++) {
            // loops over the transverse plane
            double x_local = grid_x0 + (i + 0.5)*grid_dx;
            for (int j = 0; j < ny; j++) {
                double y_local = grid_y0 + (j + 0.5)*grid_dy;
                bool intersect = check_intersect(T_cut, tau_local, x_local,
                                                 y_local, grid_dt, grid_dx,
                                                 grid_dy, cube);
                if (intersect) {
                    // cube already has the required std::array layout; pass it directly
                    cornelius_ptr->find_surface_3d(cube);
                    for (int isurf = 0; isurf < cornelius_ptr->get_number_elements();
                         isurf++) {
            double tau_center = (
                cornelius_ptr->get_centroid_element(isurf, 0)
                + tau_local - grid_dt/2.);
            double x_center = (
                cornelius_ptr->get_centroid_element(isurf, 1)
                + x_local - grid_dx/2.);
            double y_center = (
                cornelius_ptr->get_centroid_element(isurf, 2)
                + y_local - grid_dy/2.);

            double da_tau =
                cornelius_ptr->get_normal_element(isurf, 0);
            double da_x = cornelius_ptr->get_normal_element(isurf, 1);
            double da_y = cornelius_ptr->get_normal_element(isurf, 2);
                       
                        if (hydro_type == 1) {
#ifdef USE_HDF5
                            hydroinfo_ptr->getHydroinfo(
                                tau_center, x_center, y_center, fluidCellptr);
#endif
                        } else {
                            hydroinfo_MUSIC_ptr->getHydroValues(
                                x_center, y_center, 0.0, tau_center,
                                fluidCellptr);
                        }

                        output << scientific << setw(18) << setprecision(8) 
                               << tau_center << "   " << x_center << "   "
                               << y_center << "   " 
                               << da_tau << "   " << da_x << "   "
                               << da_y << "   " 
                               << fluidCellptr->temperature << "   "
                               << fluidCellptr->vx << "   " << fluidCellptr->vy 
                               << endl;

                    }
                }
            }
        }
    }
    output.close();
    
    delete fluidCellptr;
    return 0;
}
