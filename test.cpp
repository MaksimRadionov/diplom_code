#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <vector>
#include <unistd.h>
#include <stdint.h>
#include <string>
#include <cmath>
// size of double - 8 bytes always
// long - mb 4 mb 8
//fseek(File, 9 bute,SEEK_CUR(SET,END))
size_t LONG_SIZE = 4;
size_t DOUBLE_SIZE = 8;
struct InitParameters
{
    int version, wigth, height, sample_lenght;
    double horz_step, vert_step, freq,zero_pos_1{1},zero_pos_2{1};
    InitParameters();
};

InitParameters::InitParameters()
{
    fread(&version, LONG_SIZE, 1, stdin);
    fread(&wigth, LONG_SIZE, 1, stdin);
    fread(&height, LONG_SIZE, 1, stdin);
    fread(&horz_step, DOUBLE_SIZE, 1, stdin);
    fread(&vert_step, DOUBLE_SIZE, 1, stdin);
    fread(&sample_lenght,LONG_SIZE, 1, stdin);
    fread(&freq, DOUBLE_SIZE, 1, stdin);
}

double contrast(std::vector<double>& payload, double max, int left, int right)
{
    double contrast = 0;
    for(auto it = payload.begin()+left; it <= payload.begin() + right; it++)
        contrast += (*it) * (*it);
    contrast /= (right - left +1) * (right - left );//n = r-l +1 (?)
    contrast = sqrt(contrast);
    contrast = max / contrast;
    return contrast;
}

double by_3_points(std::vector<double>& y, int n)
{
    return n + (y[n-1] - y[n+1])/(2 *(y[n-1] - 2  * y[n] + y[n+1]));
}
int main(int argc, char* argv[])
{
    FILE* file = fopen("graph.txt","w");
    FILE* file2 = fopen("graph2.txt","w");
    std::vector<double> payload = {1,1,1,3,2,1,1,2,1,4,2,1};
//                                   |         |
    //    std::vector<double> payload = {1,2,3,4,5,6,7,8,9,4,2,1};
    
    int bottom_time = 4;
    double coeff_1 = 0.6;
    double coeff_2 = 1.8;
    int left = bottom_time * coeff_1 - 1;//2-1!!!!!!!!!!!!!!!
    int right = bottom_time * coeff_2 - 1;//7-1 
    
    
    double max = 0; 
    for (auto it = payload.begin(); it!=payload.end(); ++it)
    {
        fprintf(file,"%f\n",*it);
    } 
    //max_n = std::max_element(payload.begin(),payload.end()) - payload.begin();
    std::vector<double>::iterator max_n = 
            std::max_element(payload.begin() + left ,payload.begin() + right); 
    max = *max_n;
    int n = max_n - payload.begin();
    fprintf(stderr,"contrast %f\n",contrast(payload, max, left,right));// nu i dalshe
    double max_dob = by_3_points(payload,n);    
    
    fprintf(stderr, "max_n %d max %f max_d %f left %f right %f \n", n,max,max_dob, *(payload.begin() + left),*(payload.begin() + right));
//    fprintf(file2,"%d %f\n%d %f\n%d %f\n%d %f\n%d %f\n%d %f",left, 0.0, left, 0.07, n, max, right, 0.07, right, 0.0, n, 0.0);
    
    fclose(file);
    return 0;
}




