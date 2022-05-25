#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <fstream>
#include <set>
#include <iomanip>
#include <sstream>
using namespace std;

int main()
{
    set<unsigned> ns = {10};
    //unsigned n = 2;
    double b = 0.1;
    size_t L = 16;
    int lambda = 1000;
    size_t R = 1;
    for(auto n : ns) {
        ofstream fout("F:\\Ising_Kasimir\\output\\out_lambda_" + to_string(lambda) + "_R_" + to_string(R) + "_n_" + to_string(n) + ".txt");
        for(b=0.1; b<=4.1; b+=0.1) {
            fout << b;
            stringstream ss;
            ss << std::fixed << std::setprecision(1) << b;
            string pathin = "F:\\Ising_Kasimir\\observs\\obs_L_" + std::to_string(L) + "_n_"
                + std::to_string(n) + "_b_" + ss.str() + "_lambda_" + std::to_string(lambda) + "_R_" + std::to_string(R) + ".txt";
            ifstream fin(pathin);
            string s;
            for(int i = 0; i<99; ++i)
                getline(fin, s); 
            int k;
            vector<double> plaqs;
            while(fin >> k) {
                double plaq;
                for (int i = 0; i < 5; ++i) {
                    fin >> plaq;
                    plaqs.push_back(plaq);
                }
            }
            
            
            for (int i = 0; i < 5; ++i) {
                double sum = 0;
                for (k = 0; k < plaqs.size()/5.0; ++k) {
                    sum += plaqs[k*5+i];
                }
                double aver = sum / double(plaqs.size()/5.0);
                
                fout << "\t" << aver;
            }
            fout << "\n";
            fin.close();
        }
        fout.close();
    }
    return 0;
}
