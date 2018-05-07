#include <cstdlib>
#include <cstdio>
#include <vector>
#include <unistd.h>
#include <stdint.h>
#include <string>
#include <algorithm>
#include <cmath>
#include "compensation.hpp"
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

    std::vector<double> payload;
    Compensator comp(4000); 
    Filter filter(4096);
    comp.init_signal_and_spectrum();
    //comp.difraction(9);
    //comp.write_init_and_difr();
    FILE * file_nc = fopen("matrix_n_comp.txt", "w");
    FILE * file_c = fopen("matrix_comp.txt", "w");
    FILE * file_fc = fopen("matrix_fcomp.txt", "w");
    InitParameters init_parameters;
    fprintf(stderr,"wigth %d\n", init_parameters.wigth);
    fprintf(stderr,"height %d\n",init_parameters.height);
    int current_wigth =1;
    int current_height =1;
    int bottom_time = atoi(argv[1]);
    fprintf(stderr,"botom tim %d\n", bottom_time);
    double coeff_1 = 0.97;
    double coeff_2 = 1.03;
    int left = bottom_time * coeff_1;
    int right = bottom_time * coeff_2;//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!границы поиска максимума
    //int left = bottom_time -60;
    //int right = bottom_time * +60;
    double max =0;
    /*double** matrix;
    matrix = new double* [init_parameters.height];
    for(int i=0; i<init_parameters.height;i++)
            matrix[i] = new double[init_parameters.wigth];
    */
    long double sredn = 0;
    for(current_height = 1; current_height <= init_parameters.height; current_height++)
    {       
        for(current_wigth = 1; current_wigth <= init_parameters.wigth; current_wigth++)
        {
            payload.resize(4096);
            fread(&init_parameters.zero_pos_1, DOUBLE_SIZE, 1, stdin);
            if(init_parameters.zero_pos_1)
                fprintf(stderr, "ALARM %f \n", init_parameters.zero_pos_1);
            fread(&payload[0], DOUBLE_SIZE, 4096, stdin);

            filter.make_filtering(payload,1000000,6000000);
            
            std::vector<double>::iterator max_n = (std::min_element(payload.begin() + left ,payload.begin() + right));
            max = *max_n;
            int n = max_n - payload.begin();
            //fprintf(stderr,"%f \n",contrast(payload, max, left, right));
            //if(contrast(payload, max, left, right) < 10)
            //       fprintf(stderr,"contrast < 10");
            double non_comp_max =   by_3_points(payload,n);           
            auto dif_tuple = comp.get_comp_max(non_comp_max);
            fprintf(file_nc,"%f ", /*3000000.0/*/non_comp_max);
            //if (non_comp_max > 1044.0)
              //  fprintf(stderr,"h%d w %d  %f \n", current_height,current_wigth,non_comp_max);
            sredn+=non_comp_max;//+ std::get<0>(dif_tuple);
            //fprintf(stderr,"comp doba %lf\n",(double) std::get<1>(dif_tuple));
            /*if (current_height==12)
            {fprintf(stderr,"h %d w %d 3_p %f  n %d p[n-1] %f p[n]%f p[n+1] %f \n", current_height,current_wigth,non_comp_max, n, payload[n-1],payload[n],payload[n+1]);
                
            }*/
            fprintf(file_fc,"%f ", /*3000000.0/*/(non_comp_max + std::get<0>(dif_tuple)));
            fprintf(file_c,"%f ", /*3000000.0/*/(non_comp_max + std::get<1>(dif_tuple)));
             
            //fprintf(stdout,"non %f difr %f full %f \n",non_comp_max,non_comp_max + std::get<1>(dif_tuple), non_comp_max + std::get<0>(dif_tuple));
            //if (n<2970)
              //  fprintf(stderr,"h %d w %d max %d\n",current_height,current_wigth, n);

            //fprintf(stderr,"%f \n", by_3_points(payload,n));
        }
        fprintf(file_nc,"\n");
        fprintf(file_c,"\n");
        fprintf(file_fc,"\n");
    }
    fclose(file_nc);
    fclose(file_c);
    fclose(file_fc);
    sredn=sredn/(init_parameters.height*init_parameters.wigth);
   //sredn = 3000000.0/sredn; 
   fprintf(stderr,"srednyaa skorost %lf\n",(double) sredn);
    FILE * file_nc2 = fopen("matrix_n_comp.txt", "r");
    FILE * file_ncaver = fopen("matrix_comp-a.txt", "w");
    double buf;
    double max_dif = 0;
    for(current_height = 1; current_height <= init_parameters.height; current_height++)
    {       
        for(current_wigth = 1; current_wigth <= init_parameters.wigth; current_wigth++)
        {
            fscanf(file_nc2,"%lf",&buf);
            double dif = (double)sredn - buf;
            fprintf(file_ncaver,"%lf ", dif);
            if (fabs(dif) > max_dif)
                max_dif = fabs(dif);
        }
            fprintf(file_ncaver,"\n");
    }
    fclose(file_nc2);
    
   fprintf(stderr,"max dif  %lf\n",(double) max_dif);
    return 0;
}




