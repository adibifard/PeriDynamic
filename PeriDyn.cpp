#include <hpx/hpx.hpp>
#include <hpx/modules/format.hpp>
#include <hpx/debugging/print.hpp>
#include <hwloc/helper.h>
#include <hpx/include/parallel_for_loop.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/include/lcos.hpp>
#include <hpx/include/parallel_algorithm.hpp>
#include <hpx/include/async.hpp>
#include <hpx/include/util.hpp>

#include <hpx/hpx_main.hpp>

#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <math.h>
#include <string>

//Task 1: Data structures
//Implement a struct/class to store three-dimensional vectors
//Overload the operators and functions you need
//Add all the std::vector objects you need to store the simulation data
//Add all other variables you need to control the simulation
//Task 1: Data structures
//Implement a struct/class to store three-dimensional vectors
//Overload the operators and functions you need
//Add all the std::vector objects you need to store the simulation data
//Add all other variables you need to control the simulation
template<typename T>
struct vector {
    T x, y, z;

    T norm() { return std::sqrt(x * x + y * y + z * z); };

    // Vector construction
    vector(T x = 0, T y = 0, T z = 0) : x(x), y(y), z(z) {}

    // Overload + operator
    vector<T> operator+(const vector<T> rhs) {
        return vector<T>(x + rhs.x, y + rhs.y, z + rhs.z);
    }

    // Overload - operator
    vector<T> operator-(const vector<T> rhs) {
        return vector<T>(x - rhs.x, y - rhs.y, z - rhs.z);
    }

    // Overload the == operator
    bool operator==(const vector<T> rhs) const {
        return x == rhs.x
            && (y == rhs.y)
            && (z == rhs.z);
    }

    // Overload * for multiplication of vector with a scalar
    vector<T> operator*(const double rhs) {
        return vector<T>(x * rhs, y * rhs, z * rhs);
    }
    // Return size of the vector
    //std::size_t size() const { return size; }

    // Return the data of an index
    double operator[](std::size_t idx) {
        if (idx == 1) {
            return x;
        }
        else if (idx == 2)
        {
            return y;
        }
        else if (idx == 3) {
            return z;
        }
    }
    friend std::istream& operator>>(std::istream&, vector<T>&);
};

extern "C++"
std::istream & operator>>(std::istream & in, vector<double> & v) {
	return in >> v.x >> v.y >> v.z;
}

extern "C++"
std::ostream & operator<<(std::ostream & out, vector<double> & v) {
	auto precision = out.precision();
	auto width = out.width();
	out << std::setw(width) << std::setprecision(precision) << v.x << "  ";
	out << std::setw(width) << std::setprecision(precision) << v.y << "  ";
	out << std::setw(width) << std::setprecision(precision) << v.z;
	return out;
}

class PreDyn {
private:
    double delta; // Neighborhood radius (from config.dat)

    //Dimensions of the box
    double lx = 1;
    double ly = 1;
    double lz = 1;
    double h = 0.1;//mesh size

    double c;
    double sc;

    // Define the needed vectors in the code
    std::vector<double> dV; //volume of each element - To be calculated within the code

    std::vector<vector<double>> forces; // The force vector 
    std::vector<vector<double>> b; // The external force or the boundary conditions
    //std::vector<std::vector<int>> NeighbList; // List of the neighborhood

public:

    // Calculate nx,ny,nz
    size_t nx = lx / h;
    size_t ny = ly / h;
    size_t nz = lz / h;
    size_t n = nx * ny * nz; //total number of nodes

    // A method to calculate the initial position and volume of each node
    // Input > lx,ly,lz and h
    void CalcPosVol() {
        size_t ind = 0;
        vector<double> NodP;
        vector<double> B;
        double Bc = 0.1;
        for (size_t i = 0; i < nx; ++i) {
            for (size_t j = 0; j < ny; ++j) {
                for (size_t k = 0; k < nz; ++k) {
                    NodP = { h / 2 + k * h,h / 2 + j * h,h / 2 + i * h };
                    positions[ind] = positions[ind] + NodP;
                    // Assign b-values (External force-BCs)
                    if (k == 0) {  
                        B={ -Bc,0.0,0.0};
                        b[ind] = B;
                    }
                    else if (k == nx - 1) {
                        B = {Bc,0.0,0.0};
                        b[ind] = B;
                    }
                    else {
                        B = {0.0,0.0,0.0};
                        b[ind] = B;
                    }
                    ind += 1;
                }
            }
        }
    }

    //////////////////////// MATERIAL MODEL ////////////////////////////////
    // Write a function which computes the neighborhoods and store the results in a std::vector.
    std::vector<int> Neighb(size_t ni) {
        vector<double> diff;
        std::vector<int> NeighbList;
        for (size_t i = 0; i < n; ++i) {
            if (i != ni) {
                diff = positions[i] - positions[ni];
                if (diff.norm() <= delta) {
                    NeighbList.push_back(i);
                }
            }
        }
        return NeighbList;
    }

    // Write a function which computes the bond stretch
    double CompBondStretch(size_t i, size_t j) {
        vector<double> distTmp;
        vector<double> diffDispl;
        vector<double> dist_diffDispl;

        // The difference between position of the nodes
        distTmp = positions[i] - positions[j];
        // Calculate the difference in the current displacement (which is the 
        diffDispl = dist1[i] - dist1[j];

        dist_diffDispl = distTmp + diffDispl;
        // Calculate bond stretch
        double s = (dist_diffDispl.norm() - distTmp.norm()) / distTmp.norm(); // The code needs to return this vallue
        return s;
    }

    // Write a function which computes the pair-wise force
    vector<double> CompPairForce(size_t i, size_t j) {
        vector<double> PairForce;
        double s = CompBondStretch(i, j);
        double mu = CompDamage(i, j);

        vector<double> distTmp;
        vector<double> diffDispl;
        vector<double> dist_diffDispl;

        distTmp = positions[i] - positions[j];
        diffDispl = dist[i] - dist[j];
        dist_diffDispl = distTmp + diffDispl;

        PairForce = dist_diffDispl * (c * s * mu / dist_diffDispl.norm());
        return PairForce;
    }

    // Write a function which computes the damage of the node $X_i$ 
    int CompDamage(size_t i, size_t j) {
        int mu;
        double s = CompBondStretch(i, j);
        if (s < sc) {
            mu = 1;
        }
        else {
            mu = 0;
        }
        return mu;
    }

    ////////// Compute the forces for each node  ///////////////////////////////
    // A function to calculate forces for each partition
    std::vector<vector<double>> force_partition(std::vector<vector<double>> partition_data,
        std::vector<size_t> part_ind) {

        std::vector<vector<double>> ff(part_ind.size());
        for (size_t i = 0; i < part_ind.size(); ++i) {
            vector<double> sumF = { 0,0,0 };
            // Find the neighberhood of i
            std::vector<int> NeighbList;
            NeighbList = Neighb(part_ind[i]);
            // Use parallelization to claculate the pair of forces
            hpx::parallel::for_loop(hpx::parallel::execution::par, 0, NeighbList.size(),
                [&](size_t j) {
                        // Calculate the pairwise forces
                        vector<double> PairF = CompPairForce(part_ind[i], NeighbList[j]); //vector<double>
                        sumF = sumF + PairF; // Sum all the forces
                    });
            ff[i] = sumF + b[part_ind[i]];
        }
        return ff;
    }


    void computeForces() {
        size_t Ng = ceil(n/10); // It specifies the number of points calculated per node (length of each subvector)
        // Determine the number of groups for parallel calculations
        size_t const Nsv = (n - 1) / Ng + 1;

        // an index vector
        std::vector<size_t> NInd(n);
        std::generate(NInd.begin(), NInd.end(), [n = 0]() mutable { return n++; });

        std::vector<std::vector<vector<double>>> Vsub;// A 3D vector to keep partition data
        Vsub.resize(Nsv);
        
        std::vector<std::vector<size_t>> VsubIndex; // A 2D vector to trace indices of the data
        VsubIndex.resize(Nsv);

        // Partition the position data into subvectors and memorize the indices
        for (size_t i = 0; i < Nsv; ++i) {
            size_t start_ind =i*Ng;
            size_t end_itr =i*Ng+Ng;
            Vsub[i].resize(Ng);
            // handle the last subvector
            if (i * Ng + Ng > positions.size()) {
                end_itr = positions.size();
                Vsub[i].resize(positions.size() - i * Ng);
            }
            // Copy elements into the 3D subvector for Vsub and 2D vector for VsubIndex
            size_t k = 0;
            for (size_t j = start_ind; j < end_itr; j++) {
                Vsub[i][k] = positions[j];
                VsubIndex[i].push_back( NInd[j]);
                k += 1;
            }
        }

        auto fp = hpx::util::unwrapping(&PreDyn::force_partition);
        
        // Pass each subvector into the force_partition function
        // Nsv is the number of subvectors
        std::vector<hpx::lcos::future<std::vector<vector<double>>>> F;
        F.resize(Nsv);
       for (size_t i = 0; i < Nsv; ++i) {
          F[i]=hpx::async(&PreDyn::force_partition,this, std::move(Vsub[i]), std::move(VsubIndex[i]));
        }
       // Map data back into the forces vector
       std::vector<vector<double>> fAll;
       auto f=hpx::when_all(F).then([&fAll](auto&& f) {
           auto futures = f.get();
           size_t k = 0;
           for (size_t i = 0; i < futures.size(); i++) {
               auto h = futures[i].get();
               for (size_t j = 0; j < h.size(); j++) {
                   fAll.push_back(h[j]);
                   k += 1;
               }
           }
           });
       f.get();
       // Collect all forces into the forces vector
       for (size_t i = 0; i < fAll.size(); i++) {
           forces[i] = fAll[i];
       }
    }

    // A function to update the displacement
    void UpdateDisp() {
        // Update the displacement vector using the hpx::parallel::for_loop
        hpx::parallel::for_loop(hpx::parallel::execution::par, 0, forces.size(),
            [&](size_t i) {dist1[i] = dist[i] * 2 - disOld[i] + forces[i] * pow(timeStepSize, 2);});
        std::copy(dist.begin(), dist.end(), disOld.begin());
        std::copy(dist1.begin(), dist1.end(), dist.begin());
    }

    // A function to update the position
    void UpdatePosition() {
        // Update the displacement vector using the hpx::parallel::for_loop
        hpx::parallel::for_loop(hpx::parallel::execution::par, 0, positions.size(),
            [&](size_t i) {positions[i] = positions[i] + dist1[i]; });
    }

public:
    std::vector<vector<double>> positions;
    std::vector<vector<double>> dist1; //future (t+1) displacement
    std::vector<vector<double>> dist; //Current (t) displacement
    std::vector<vector<double>> disOld;//previous (t-1) displacement


    double K; //Material property (from config.dat)
    double K_ic; //Material property (from config.dat)
    double ro; // Mass Density (from config.dat)
    size_t tf; //final time (from config.dat)
    size_t timeStepSize; //Time-Step Size (from config.dat)
    size_t timeSteps = 3;
    //Constructor
    PreDyn(std::string& fileName) {
        std::ifstream ifs(fileName);
        if (!ifs.is_open()) {
            throw std::runtime_error("Could not open " + fileName + "!");
        }
        /*Read the material properties $(K,K_{Ic},\varrho)$, the final time T,
        the time step with $t_n$, and the horizon $\delta$ */
        ifs >> K >> K_ic >> tf >> timeStepSize ;
        delta = 3*h;
        // Calculate c (Stiffness) and Sc (Critical bond strength)
        c = (18 * K) / (3.14 * delta);
        sc = (5 / 12) * sqrt(K_ic / (pow(K, 2) * delta));
        sc = 0.04;
        c = 100;
        std::cout << c;
        dV.resize(n);
        std::fill(dV.begin(), dV.end(), double());

        positions.resize(n);
        std::fill(positions.begin(), positions.end(), vector<double>());

        dist.resize(n);
        std::fill(dist.begin(), dist.end(), vector<double>(0.0));

        disOld.resize(n);
        std::fill(disOld.begin(), disOld.end(), vector<double>(0.0));

        dist1.resize(n);
        std::fill(dist1.begin(), dist1.end(), vector<double>());

        forces.resize(n);
        std::fill(forces.begin(), forces.end(), vector<double>());

        b.resize(n);
        std::fill(b.begin(), b.end(), vector<double>());

        // Initialize the position vector
        CalcPosVol();
        // Print the data to the Console
        std::cout << "Contents of " << fileName << std::endl;
        std::cout << K << ' ' << K_ic << ' ' << tf << ' ' << timeStepSize << ' ' << delta << std::endl;
        /*for (int i = 0; i < n; ++i) {
            std::cout << positions[i] << std::endl;
        }*/
        std::cout << std::endl << "Data   :      x          y          z    |     vx         vy         vz" << std::endl;

    }

    // A method to save the output file
    void saveFile(std::vector<vector<double>> Data,std::string name) {
        std::ofstream myfile;
        std::string ext = ".csv";
        std::string fullName = name + ext;
        myfile.open(fullName);

        for (size_t j = 0; j < Data.size(); j++) {
                myfile << Data[j].x <<","<< Data[j].y<<","<< Data[j].z<< std::endl;
            }
        myfile.close();
    }
    friend std::ostream& operator<<(std::ostream&, PreDyn&);
};

extern "C++"
std::ostream & operator<<(std::ostream & out, PreDyn & pd) {
    for (int i = 0; i < pd.n; ++i) {
        out << "Element " << i + 1 << " : ";
        out << std::setprecision(6) << std::setw(9) << pd.dist1[i];
        out << std::endl;
    }
    return out;
}


int main()
{
    std::string Input= "config.dat";
    std::string fnameDist = "Displacement_";

    // instantiate an instance of the PreDyn class
    PreDyn pd(Input);

    for (size_t t = 0; t < pd.timeSteps; t++) {
        std::cout << std::endl << "Cycle " << t << std::endl;
        // Determine the force-field
        pd.computeForces();
        // Update the displacement 
        pd.UpdateDisp();
        // Update the positions
        pd.UpdatePosition();
        // Write dist(x,y,z) to the output file
        std::string Tstep = std::to_string(t);
        std::string fileName = fnameDist + Tstep;
        pd.saveFile(pd.dist1, fileName);

        std::cout << pd;
    }
    return EXIT_SUCCESS;

}
