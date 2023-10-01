#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main() {
    // Simulation Parameters
    cout << "Start \n";
    clock_t start, end;

    const int Nx = 51;
    const int Ny = 51;
    const int Re = 100;
    const double tol = 1e-12;
    const double alpha = 0.8;
    const double Lx = 1;
    const double Ly = 1;
    string name = "Steady Lid Driven Flow - Re "+to_string(Re)+" - "+to_string(Nx)+"x"+to_string(Ny);
    cout << name;

    // Derived Parameters
    double dx = Lx/(Nx-1);
    double dy = Ly/(Ny-1);

    // Allocate Space for Variables
    double u_saved[Ny][Nx] = {0};
    double v_saved[Ny][Nx] = {0};
    double p_saved[Ny][Nx] = {0};

    // Staggered Variables
    double u[Ny+1][Nx] = {0};
    double u_star[Ny+1][Nx] = {0};
    double u_m[Ny+1][Nx] = {0};
    for(int i = 0; i < Nx; ++i) {
        u[0][i] = 2;
        u_star[0][i] = 2;
        u_m[0][i] = 2;
    }
    double d_e[Ny+1][Nx] = {0};

    double v[Ny][Nx+1] = {0};
    double v_star[Ny][Nx+1] = {0};
    double v_m[Ny][Nx+1] = {0};
    double d_n[Ny][Nx+1] = {0};

    double p[Ny+1][Nx+1] = {0};
    double p_m[Ny+1][Nx+1] = {0};
    double pc[Ny+1][Nx+1] = {0};
    double b[Ny+1][Nx+1] = {0};

    // Solve Governing Equations
    int iter = 0;
    double error = 1;

    double u_e, u_w, v_n, v_s, A_E, A_W, A_N, A_S, A_P;
    start = clock();
    while(error > tol) {
        // X-Momentum Equation
        for(int i = 1; i < Ny; ++i) {
            for(int j = 1; j < Nx-1; ++j) {
                cout << i;
                return 0;
                u_e = (u[i][j] + u[i][j+1])/2;
                u_w = (u[i][j] + u[i][j-1])/2;
                v_n = (v[i-1][j] + v[i-1][j+1])/2;
                v_s = (v[i][j] + v[i][j+1])/2;

                A_E = u_e*dy/2-dy/(Re*dx);
                A_W = -u_w*dy/2-dy/(Re*dx);
                A_N = v_n*dx/2 - dx/(Re*dy);
                A_S = -v_s*dx/2 - dx/(Re*dy);

                A_P = (u_e*dy-u_w*dy+v_n*dx-v_s*dx)/2+2*dy/(Re*dx)+2*dx/(Re*dy);
                d_e[i][j] = -dy/A_P;
                u_star[i][j] = -(A_E*u[i][j+1]+A_W*u[i][j-1]+A_N*u[i-1][j]+A_S*u[i+1][j])/A_P - d_e[i][j] * (p[i][j] - p[i][j+1]);
            }
        }

        // X Momentum Boundary Conditions
        for(int i = 0; i < Nx;++i) { // Top and Bottom BC
            u_star[0][i] = 2 - u_star[1][i];
            u_star[Ny][i] = -u_star[Ny-1][i];
        }
        for(int i = 1; i < Ny; ++i) {// Left and Right BC
            u_star[i][0] = 0;
            u_star[i] [Nx-1] = 0;
        }

        // Y Momentum
        for(int i = 1; i < Ny-1;++i) {
            for(int j = 1; j < Nx;++j) {
                u_w = (u[i][j-1] + u[i+1][j-1])/2;
                u_e = (u[i][j] + u[i+1][j])/2;
                v_n = (v[i][j] + v[i-1][j])/2;
                v_s = (v[i][j] + v[i+1][j])/2;

                A_E = u_e*dy/2 - dy/(Re*dx);
                A_W = -u_w*dy/2 - dy/(Re*dx);
                A_N = v_n*dx/2 - dx/(Re*dy);
                A_S = -v_s*dx/2 - dx/(Re*dy);

                A_P = (u_e*dy - u_w*dy + v_n*dx - v_s*dx)/2 + 2*dy/(Re*dx) + 2*dx/(Re*dy);
                d_n[i][j] = -dx/A_P;

                v_star[i][j] = -(A_E*v[i][j+1] + A_W*v[i][j-1] + A_N*v[i-1][j] + A_S*v[i+1][j])/A_P - d_n[i][j] * (p[i+1][j] - p[i][j]);
            }
        }

        // Y Momentum Boundary Conditions
        for(int i = 0;i < Ny;++i) {// Left and Right BC
            v_star[i][1] = -v_star[i][2];
            v_star[i][Nx] = -v_star[i][Nx-1];
        }
        for(int i = 0;i < Nx;++i) {// Top and Bottom BC
            v_star[1][i] = 0;
            v_star[Ny-1][i] = 0;
        }

        // Pressure Correction
        for(int i = 1; i < Ny+1; ++i) {
            for(int j = 1; j < Nx+1; ++j) {
                A_E = -d_e[i][j]*dy;
                A_W = -d_e[i][j-1]*dy;
                A_S = -d_n[i][j]*dx;
                A_N = -d_n[i-1][j]*dx;
                A_P = A_E+A_W+A_S+A_N;
                
                b[i][j] = -(u_star[i][j]-u_star[i][j-1])*dy+(v_star[i][j]-v_star[i-1][j])*dx;
                pc[i][j] = b[i][j]/A_P;
            }
        }

        // Apply Pressure Correction
        for(int i = 1;i<Ny;++i) {
            for(int j = 1;j<Nx;++j) {
                p_m[i][j] = p[i][j] + alpha*pc[i][j];
            }
        }

        // Pressure Boundary
        for(int i = 0;i <Nx+1;++i){// Top and Bottom BC
            p_m[0][i] = p_m[1][i];
            p_m[Ny][i] = p_m[Ny-1][i];
        }
        for(int i = 0; i < Ny+1;++i){// Left and Right BC
            p_m[i][0] = p_m[i][1];
            p_m[i][Nx] = p_m[i][Nx-1];
        }

        // Correct Velocities
        // X Correction
        for(int i = 1;i<Ny;++i){
            for(int j = 1;j<Nx-1;++j) {
                u_m[i][j] = u_star[i][j]+alpha*d_e[i][j]*(pc[i][j+1]-pc[i][j]);
            }
        }

        // X Boundary
        for(int i = 0; i < Nx;++i) { // Top and Bottom BC
            u_m[0][i] = 2 - u_m[1][i];
            u_m[Ny][i] = -u_m[Ny-1][i];
        }
        for(int i = 1; i < Ny; ++i) {// Left and Right BC
            u_m[i][0] = 0;
            u_m[i][Nx-1] = 0;
        }

        // Y Correction
        for(int i = 1;i<Ny-1;++i) {
            for(int j = 1;j<Nx;++j) {
                v_m[i][j] = v_star[i][j] + alpha*d_n[i][j]*(pc[i][j]-pc[i+1][j]);
            }
        }

        // Y Boundary
        for(int i = 0;i < Ny;++i) {// Left and Right BC
            v_m[i][0] = -v_m[i][1];
            v_m[i][Nx] = -v_m[i][Nx-1];
        }
        for(int i = 0;i < Nx;++i) {// Top and Bottom BC
            v_m[0][i] = 0;
            v_m[Ny-1][i] = 0;
        }

        // Update Error
        error = 0;
        for(int i = 1;i<Ny;++i){
            for(int j = 1;j<Nx;++j) {
                error += abs(b[i][j]);
            }
        }
        double du = 0;
        
        // Update User
        if(iter%1000 == 0) {
            cout << "Iteration: ";
            cout << iter;
            cout << " Error: ";
            cout << error << "\n";

        }
        
        // Update Variables
        for(int i = 0;i<Ny+1;++i) {
            for(int j = 0;j<Nx;++j) {
                u[i][j] = u_m[i][j];
            }
        }
        for(int i = 0;i<Ny;++i) {
            for(int j = 0;j<Nx+1;++j) {
                v[i][j] = v_m[i][j];
            }
        }
        for(int i = 0;i<Ny+1;++i) {
            for(int j = 0;j<Nx+1;++j) {
                p[i][j] = p_m[i][j];
            }
        }
        ++iter;
    }
    end = clock();
    // Evaluate Colocated Results
    double xset[Ny][Nx] = {0};
    double yset[Ny][Nx] = {0};
    double x = 0;
    double y = Ly;
    for(int i = 0;i<Ny;++i) {
        x = 0;
        for(int j = 0;j<Nx;++j) {
            u_saved[i][j] = (u_m[i][j]+u_m[i+1][j])/2;
            v_saved[i][j] = (v_m[i][j]+v_m[i][j+1])/2;
            p_saved[i][j] = (p_m[i][j]+p_m[i][j+1]+p_m[i+1][j]+p_m[i+1][j+1])/4;

            xset[i][j] = x;
            yset[i][j] = y;

            x += dx;
        }
        y -= dy;
    }
    double dt = double(end-start);
    dt /= double(CLOCKS_PER_SEC);
    cout << "\nTime Taken: " << dt;
    // Write Results to CSV file
    // Results organized as follows:
    // x,y,u,v,p
    // First row tells us dx,dy,tol,tol,tol
    ofstream file;
    file.open(name+".csv");
    file << dx << ",";
    file << dy << ",";
    file << tol << ",";
    file << tol << ",";
    file << tol << "\n";
    for(int i = 0; i <Ny;++i) {
        for(int j = 0;j<Nx;++j) {
            file << xset[i][j] << ",";
            file << yset[i][j] << ",";
            file << u_saved[i][j] << ",";
            file << v_saved[i][j] << ",";
            file << p_saved[i][j] << "\n";
        }
    }
    file.close();
    return 0;
}