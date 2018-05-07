#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <vector>
#include <unistd.h>
#include <stdint.h>
#include <string>
#include <cmath>
#include "filter.hpp"
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
    std::string file_name = argv[1];

    std::vector<double> payload;
    
    FILE * file = fopen("graph.txt", "w");
    FILE * file2 = fopen("graph2.txt", "w");
    InitParameters init_parameters;
    Filter filter(4096); 
    fprintf(stderr,"height %d\n",init_parameters.height);
    fprintf(stderr,"wigth %d\n", init_parameters.wigth);
    int current_wigth ;
    int current_height ;
    current_wigth = atoi(argv[2]);
    current_height = atoi(argv[1]);
    int bottom_time = atoi(argv[3]);
    double coeff_1 = 0.95;
    double coeff_2 = 1.05;
    int left = bottom_time * coeff_1;
    int right = bottom_time * coeff_2;//!!!!!!!!!!!!!!!!!!!!!границы поиска максимума
    //int left = bottom_time - 50;
    //int right = bottom_time + 50;
    
    if((current_wigth > init_parameters.wigth || current_wigth < 1)||
           (current_height > init_parameters.height || current_height < 1))
    {
        fprintf(stderr, "Wrong parameters\n");
        return 1;
    }
    fprintf(stderr,"current_height %d\n",current_height);
    fprintf(stderr,"current_wigth %d\n",current_wigth);
    
    
    fseek(stdin, ((current_height - 1) * init_parameters.wigth + (current_wigth-1)) * 
                    (init_parameters.sample_lenght + 1) * DOUBLE_SIZE, SEEK_CUR);
    payload.resize(4096);//1111111111111111111111111111111!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    fread(&init_parameters.zero_pos_1, DOUBLE_SIZE, 1, stdin);
    fread(&payload[0], DOUBLE_SIZE, 4096, stdin);
     
         filter.make_filtering(payload,1000000,6000000);

    fread(&init_parameters.zero_pos_2, DOUBLE_SIZE, 1, stdin);
    double max = 0; 
    for (auto it = payload.begin(); it!=payload.end(); ++it)
    {
        fprintf(file,"%f\n",*it);
    } 
    //max_n = std::max_element(payload.begin(),payload.end()) - payload.begin();
    std::vector<double>::iterator max_n = (std::min_element(payload.begin() + left ,payload.begin() + right));
    max = *max_n;
    int n = max_n - payload.begin();
    fprintf(stderr,"contrast %f\n",contrast(payload, max, left,right));// nu i dalshe
    fprintf(stderr, "max_n %d max %f left %f right %f \n", n,max, *(payload.begin() + left),*(payload.begin() + right));
    fprintf(file2,"%d %f\n%d %f\n%d %f\n%d %f\n%d %f\n%d %f\n%d %f",left, 0.0, left, 0.07, n, 0.06, right, 0.07, right, 0.0, n, 0.0,n, 0.06);
    fclose(file);
    return 0;
}




