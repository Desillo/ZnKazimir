#include <iostream>
#include <string>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>
#include <random>
#include "random.h"
#include <fstream>
#include <set>
#include <iomanip>
#include <sstream>

using namespace std;

#define DIM 3

class GaugeZn {
public:

    GaugeZn(size_t size, unsigned n, double beta, int lambda, size_t R, 
            bool hotstart = true) 
        : size(size), n(n), Beta(beta), Lambda(lambda), Distance(R), links(V*DIM, 0), probs(n), vcos(n), vsin(n) {
        rnd_init();
        for (unsigned i = 0; i < n; ++i)
        {
            vcos[i] = sin(ang*i);
            vsin[i] = cos(ang*i);
        }
        if(hotstart) {
            for(size_t i = 0; i < links.size(); ++i) {
                links[i] = floor(rnd_get()*double(n));
            }
        }
    }

    ~GaugeZn() {
        cout << "Computation is completed" << endl;
        rnd_free();
    }

    int& at(size_t x, size_t y, size_t z, size_t link) {
        return links[x + y*size + z*size*size + link*V]; 
    }

    int at(size_t x, size_t y, size_t z, size_t link) const {
        return links[x + y*size + z*size*size + link*V]; 
    }

    int& at(size_t x[], size_t link) {
        return links[x[0] + x[1]*size + x[2]*size*size + link*V]; 
    }

    int at(size_t x[], size_t link) const {
        return links[x[0] + x[1]*size + x[2]*size*size + link*V]; 
    }

    void HeatBath() {
        size_t x[DIM];
        size_t i,j;
        double s;
          for (x[0]=0; x[0]<size; x[0]++)
            for (x[1]=0; x[1]<size; x[1]++)
              for (x[2]=0; x[2]<size; x[2]++)
                  for (i=0; i<DIM; i++) {
                    double re_sum = 0;
                    double im_sum = 0;
                    for (j=0;j<DIM;j++){
                      if (i!=j){
                        int staple1, staple2;
                        movedown(x,j);
                        staple1 = at(x, j) - at(x, i);

                        moveup(x,i);
                        staple1 -= at(x, j);  
                        moveup(x,j);
                        // plaquette 1456 
                        staple2 = at(x, j);
                        moveup(x,j);
                        movedown(x,i);
                        staple2 -= at(x, i);
                        movedown(x,j);
                        staple2 -= at(x, j);
                        double beta = Beta;
                        if (i > 0 && j > 0 && (x[0] == 0 || x[0] == Distance))
                            beta *= Lambda;
                        
                        //staple1 = (n+sgn(staple1)*(staple1%n))%n;
                        //staple2 = (n+sgn(staple2)*(staple2%n))%n;
                        //re_sum += beta*(vcos[staple1] + vcos[staple2]);
                        //im_sum += beta*(vsin[staple1] + vsin[staple2]);
                        re_sum += beta*(cos(ang * staple1) + cos(ang * staple2));
                        im_sum += beta*(sin(ang * staple1) + sin(ang * staple2));
                      }
                    }
                    s = 0.0;
                    
                    for(unsigned k = 0; k<n; ++k) {
                        //double S = vcos[k]*re_sum - vsin[k]*im_sum;
                        double S = cos(ang*k) * re_sum - sin(ang*k) * im_sum;
                        //double S = Staples(x, i, th) - Staples(x, i, k); 
                        double p = exp(S);
                        probs[k] = p;
                        s += p;
                    }

                    double norm = double(1.0) / s;
                    double r = rnd_get();
                    for(unsigned k = 0; k < n; ++k) {
                        probs[k] *= norm;
                        if(k) probs[k] += probs[k-1];
                        if(r < probs[k]) {
                            at(x, i) = k;
                            break; 
                        }
                    }

                  }
    
    }

    double Action() const {
        size_t x[DIM];
        size_t i,j;
        double sum = 0.0;
          for (x[0]=0; x[0]<size; x[0]++)
            for (x[1]=0; x[1]<size; x[1]++)
              for (x[2]=0; x[2]<size; x[2]++)
                  for (i=0; i<DIM; i++) {
                     // cout << x[0] << " "  << x[1] << " " << x[2] << " " << i << " - " << at(x,i) << endl;
                    for (j=i+1;j<DIM;j++){
                        int staple;
                        staple = at(x, i) - at(x, j);
                        moveup(x, i);
                        staple += at(x, j);
                        moveup(x, j);
                        movedown(x, i);
                        staple -= at(x, i);
                        movedown(x, j);
                        double beta = Beta;
                        
                        //staple = (n + sgn(staple) * (staple % n)) % n;
                        //sum += vcos[staple];
                        sum += cos(ang * staple);
                    }
                  }
        return sum; 
    }
 //Compute Energy
    double Energy() const {

        return Action(); //- Action(Lambda);
    }
//Print lattice
    void PrintGzn() {
        size_t x[DIM];
        size_t i;
          for (x[0]=0; x[0]<size; x[0]++)
            for (x[1]=0; x[1]<size; x[1]++)
              for (x[2]=0; x[2]<size; x[2]++) 
                  for(i = 0; i < DIM; ++i) {
                      cout << x[0] << " "  << x[1] << " " << x[2] << " " << i << " - " << at(x,i) << endl;
                  }
                         
    }
//Anverage Plaquette
    double AverPlaquette() {
        return Energy()*on_V;
    }

    double AvPlaqOnPlates() {
        size_t x[DIM];
        size_t i, j;
        double sum = 0.0;
        for (x[0] = 1; x[0] == 0 || x[0] == Distance; x[0]+=Distance)
            for (x[1] = 0; x[1] < size; x[1]++)
                for (x[2] = 0; x[2] < size; x[2]++) {
                    int staple;
                    staple = at(x, 1) - at(x, 2);
                    moveup(x, 1);
                    staple += at(x, 2);
                    moveup(x, 2);
                    movedown(x, 1);
                    staple -= at(x, 1);
                    movedown(x, 2);
                    double beta = Beta;

                    //staple = (n + sgn(staple) * (staple % n)) % n;
                    //sum += vcos[staple];
                    sum += cos(ang * staple);
                }
        return sum/(2.0*size*size);
    }

    double AvPlaqBetweenPlates() {
        size_t x[DIM];
        size_t i, j;
        double sum = 0.0;
        for(x[0] = 0; x[0]<Distance; x[0]++)
            for (x[1] = 0; x[1] < size; x[1]++)
                for (x[2] = 0; x[2] < size; x[2]++)
                    for (i = 0; i < DIM; i++) {
                        // cout << x[0] << " "  << x[1] << " " << x[2] << " " << i << " - " << at(x,i) << endl;
                        if ((x[0] == 0) && (i > 0)) break;
                        for (j = i + 1; j < DIM; j++) {
                            int staple;
                            staple = at(x, i) - at(x, j);
                            moveup(x, i);
                            staple += at(x, j);
                            moveup(x, j);
                            movedown(x, i);
                            staple -= at(x, i);
                            movedown(x, j);
                            double beta = Beta;

                            //staple = (n + sgn(staple) * (staple % n)) % n;
                            //sum += vcos[staple];
                            sum += cos(ang * staple);
                        }
                    }
        return sum / ((3.0*Distance-1) * size * size);
    }

    double AvPlaqOutOfPlates() const {
        size_t x[DIM];
        size_t i, j;
        double sum = 0.0;
        for (x[0] = 0; x[0] < size; x[0]++) {
            if (x[0] == 0) {
                x[0] += Distance-1;

            }
            else {
                for (x[1] = 0; x[1] < size; x[1]++)
                    for (x[2] = 0; x[2] < size; x[2]++)
                        for (i = 0; i < DIM; i++) {
                            // cout << x[0] << " "  << x[1] << " " << x[2] << " " << i << " - " << at(x,i) << endl;
                            if ((x[0] == Distance) && (i > 0)) break;
                            for (j = i + 1; j < DIM; j++) {
                                int staple;
                                staple = at(x, i) - at(x, j);
                                moveup(x, i);
                                staple += at(x, j);
                                moveup(x, j);
                                movedown(x, i);
                                staple -= at(x, i);
                                movedown(x, j);
                                double beta = Beta;

                                //staple = (n + sgn(staple) * (staple % n)) % n;
                                //sum += vcos[staple];
                                sum += cos(ang * staple);
                            }
                        }
            }
        }
        return sum/((double)size * size * (3.0 * (size-Distance) - 1));
    }

    double CasimirE() {
        size_t x[DIM];
        size_t i, j;
        double sum = 0.0;
        for (x[0] = 0; x[0] < size; x[0]++) {
            if (x[0] == 0 || x[0] == Distance) continue;
            for (x[1] = 0; x[1] < size; x[1]++)
                for (x[2] = 0; x[2] < size; x[2]++) {
                    int staple;
                    staple = at(x, 1) - at(x, 2);
                    moveup(x, 1);
                    staple += at(x, 2);
                    moveup(x, 2);
                    movedown(x, 1);
                    staple -= at(x, 1);
                    movedown(x, 2);
                    double beta = Beta;

                    //staple = (n + sgn(staple) * (staple % n)) % n;
                    //sum += vcos[staple];
                    sum += cos(ang * staple);
                }
        }
        return sum / ((double)size * size * (size-2));
    }

    double PlaqOut() {
        double sum = Plaqx0(Distance, 0, 1) + Plaqx0(Distance, 0, 2);
        for(size_t x0 = Distance + 1; x0 < size; x0++)
            sum += Plaqx0(x0, 0, 1) + Plaqx0(x0, 0, 2) + Plaqx0(x0, 1, 2);
        return sum/ ((double)size * size * (3.0 * (size - Distance) - 1));
    }

    double PlaqOn() {
        double sum = Plaqx0(0, 1, 2) + Plaqx0(Distance, 1, 2);
        return sum / (2.0 * size * size);
    }

    double PlaqBetween() {
        double sum = Plaqx0(0, 0, 1) + Plaqx0(0, 0, 2);
        for (size_t x0 = 1; x0 < Distance; x0++)
            sum += Plaqx0(x0, 0, 1) + Plaqx0(x0, 0, 2) + Plaqx0(x0, 1, 2);
        return sum / ((3.0 * Distance - 1) * size * size);
    }

    double CasE() {
        double sum = 0.0;
        for (size_t x0 = 0; x0 < size; x0++)
            sum += Plaqx0(x0, 1, 2);
        return sum/((double)size*size*size);
    }

    double Plaqx0(size_t x0, size_t i, size_t j) {
        double sum = 0.0;
        size_t x[DIM];
        x[0] = x0;
        for(x[1] = 0; x[1] < size; ++x[1])
            for (x[2] = 0; x[2] < size; ++x[2]) {
                sum += Plaq(x, i, j);
            }
        return sum;
    }

    double Plaq(size_t x[], size_t i, size_t j) {
        int staple;
        staple = at(x, i) - at(x, j);
        moveup(x, i);
        staple += at(x, j);
        moveup(x, j);
        movedown(x, i);
        staple -= at(x, i);
        movedown(x, j);
        return cos(ang * staple);
    }

private:
    const size_t size;
    const unsigned n;
    const size_t V = size*size*size;
    const double on_V = double(1) / (double(V)*double(DIM));
    const double ang = double(2.0) * double(M_PI) / double(n);

    double Beta;
    int Lambda;
    size_t Distance;

    vector<int> links;
    vector<double> probs;

    vector<double> vcos;
    vector<double> vsin;

    void moveup(size_t x[], size_t d) const {
      if (x[d]==size-1) x[d]=0;
      else x[d]++;
    }

    void movedown(size_t x[], size_t d) const {
      if (x[d]==0) x[d] = size-1;
      else x[d]--;
    }

    int sgn(int val) const {
        return (val > 0) - (val < 0);
    }
/*
    double Staples(size_t x[], size_t i, int theta) {
        size_t j;
        double sum = 0.0;
        for(j = 0; j<DIM; ++j) {
            if(j==i) continue; 
            int staple1, staple2;
            movedown(x,j);
            staple1 = at(x, j) - at(x, i);

            moveup(x,i);
            staple1 -= at(x, j);  
            moveup(x,j);
            // plaquette 1456 
            staple2 = at(x, j);
            moveup(x,j);
            movedown(x,i);
            staple2 -= at(x, i);
            movedown(x,j);
            staple2 -= at(x, j);
            double B = ((x[0] == 1 || x[0] == Distance+1) && (i && j))
                ? Beta*Lambda : Beta;
            sum += B*(cos(ang*(staple1+theta)) + cos(ang*(staple2+theta)));
            
        }        
        return sum;
    }*/
};

ostream& operator<<(ostream& os, const GaugeZn& gzn) {
    os << gzn.Energy();    
    return os;
}

int main(int argc, char** argv) {
    set<unsigned> ns = {10};
    //unsigned n = 2;
    double b = 0.1;
    size_t L = 16; 
    int lambda = 1000;
    size_t R = 1;
    bool cont = false;
    unsigned MC_steps = 1000;
    for(auto n : ns) {
        for(b = 0.1; b<=4.1; b+=0.1) {
            stringstream ss;
            ss << std::fixed << std::setprecision(1) << b;
            string path = "F:\\Ising_Kasimir\\observs\\obs_L_"  + std::to_string(L) + "_n_"
                    + std::to_string(n) + "_b_" + ss.str() + "_lambda_" + std::to_string(lambda) + "_R_" + std::to_string(R) + ".txt";
            ofstream fout(path);
            fout << "#Average Plaquette and conf\n";
            fout << "#conf\tPlaq\tPlaqOn\tPlaqBetween\tPlaqOut\tCasimirE\n";
            GaugeZn gzn(L, n, b, lambda, R);
            for(unsigned k = 0; k<MC_steps; ++k) {
                gzn.HeatBath();
                double plaq = gzn.AverPlaquette();
                cout << "n = " << n << " b = " << b << " k = " << k << ": "<< plaq << endl;
                fout << k << '\t' << plaq << '\t' << gzn.PlaqOn() << '\t' << gzn.PlaqBetween() << '\t' << gzn.PlaqOut()
                    << '\t' << gzn.CasE() << '\n';
                
            }
            fout.close();
        }
    }
    return 0;
}

/*char betaname[100];
char confname[10];
strcat(betaname, "F:\\Ising_Kasimir\\observs\\obs_L_");
sprintf(confname, "%d", n);
strcat(betaname, confname);
strcat(betaname, "_n_");
sprintf(confname, "%d", n);
strcat(betaname, confname);
strcat(betaname, "_b_");
sprintf(confname, "%.1f", b);
strcat(betaname, confname);
strcat(betaname, ".txt");*/
