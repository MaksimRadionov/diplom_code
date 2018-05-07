#include "compensation.hpp"
#include <cmath>
#include <complex>
#include <algorithm>
/*
fftwf_complex* in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N);
fftwf_comple* out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * N);
*/
double cmp2 (double a,double c ,double f, double z, double f_d )//на вход частота в герцах ,f_d / f, для воввтановления передаем f_d
{
    return f_d/f;
}
double cmp (double a,double c ,double f, double z )//на вход частота в герцах ,f_d / f
{
    return c*1000*z/(pow(a,2)*M_PI*f);
}
double SQR(double a,double c ,double f, double z )
{

    return sqrt(1+pow(cmp(a,c,f,z),2));
    //return sqrt(1+pow(c*1000*z/(pow(a,2)*M_PI*f),2));
}
//компенсация с учетом размера:
//с учетом размера приемника:
complex Compensator::S(complex P0,double f0, double a,double c , double z)//на вход частота в мегагерцах
{

    if (f0>F_N) {f0 = f0 - 2*F_N;}
    double B = pow(r0,2)/(a*a*(1+pow(cmp(a,c,f0,z),2)));
    double ex1 = exp(-B);
    complex q=complex(M_PI * pow(a,2) *  (1.0 - ex1* cos(cmp(a,c,f0,z)*B)), M_PI * pow(a,2) * ex1* sin(cmp(a,c,f0,z)*B));
    return P0*q;
}
//просто дифракция:
complex Compensator::S0(complex P0, double f0, double a,double c , double z)//на вход частота в герцах
{
    if (f0>F_N) {f0 = f0 - 2*F_N;}//когда этого нет а есть f-F_N в теле то работает
    double b = 1 + pow(cmp(a,c,f0,z),2);
    
    double ksi = cmp(a,c,f0,z);//знак зависит от того че за фурье
    //double B = pow(0,2)/(a*a*b);// r0 тут стоял размерность!!
    //double ex1 = exp(-B);
    //complex q=complex(1.0 / b * (cos(ksi*B)+ksi*sin(ksi*B)), -1.0 / b * (ksi*cos(ksi*B) - sin(ksi*B)));
    complex q = complex(1.0 / b , 1.0 / b * ksi);//получается что тут  1 - i*f_d/f

    return P0*q;
}


double by_3_points(double y_m1, double y_0, double y_p1, int n)
{
    return n + (y_m1 - y_p1)/(2 *(y_m1 - 2  * y_0 + y_p1));
}

std::tuple<double,double> Compensator::get_comp_max(double non_comp_max)//тут не в микро секундах а в десятках нано надо перевести
{
   // double length = non_comp_max / 100000 * c ;   
   double length = 30 ;   
            //fprintf(stderr, "len mm   %f\n", length);
    difraction(length); 
    double new_max1 =0;
    int max_index1=0;
    double new_max2 =0;
    int max_index2=0;
    for(int i =0 ; i<Num; i++)
    {
        auto one_complex1 = P3[i];
        auto one_complex2 = P2[i];
        if (one_complex1.re>new_max1)
        {
            new_max1=one_complex1.re;
            max_index1 = i;
        }
        if (one_complex2.re>new_max2)
        {
            new_max2=one_complex2.re;
            max_index2 = i;
        }
    }
   double new_time1 = by_3_points(P3[max_index1-1].re, P3[max_index1].re, P3[max_index1+1].re, max_index1); 
   double new_time2 = by_3_points(P2[max_index2-1].re, P2[max_index2].re, P2[max_index2+1].re, max_index2); 
//    fprintf(stderr, "new_max   %f\n", new_time);
return std::make_tuple((2000 - new_time1),(2000 - new_time2));//первый полная компенсация второй только дифракция
}

void Compensator::difraction(double length)// difraction of reference signal
{
    for (int i=1; i<Num;i++)//домножаем спектр нулевую частоту пропустил
    {
        //P2[i]=S0(P_sp[i],f_[i]);// only difraction
        P2[i]=S0(P_sp[i],f_[i],a,c,length);
        P3[i]=S(P_sp[i],f_[i],a,c,length);
        //fprintf(stderr, "%d %f %f\n",i,P_sp[i].re, P_sp[i].im);
        in[i][0]=P2[i].re;
        in[i][1]=P2[i].im;
        in2[i][0]=P3[i].re;
        in2[i][1]=P3[i].im;
    }
    fftwf_plan plan1;
    plan1 = fftwf_plan_dft_1d(Num, in, out,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftwf_execute(plan1);
    fftwf_destroy_plan(plan1);
    fftwf_plan plan2;
    plan2 = fftwf_plan_dft_1d(Num, in2, out2,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftwf_execute(plan2);
    fftwf_destroy_plan(plan2);

    for(int i=0;i<Num;i++)
    {
        P2[i].re=out[i][0]/(Num);
        P2[i].im=out[i][1]/(Num);
      //  fprintf(stderr, "%d %f\n",i,out2[i][0]);
        P3[i].re=out2[i][0]/(Num);
        P3[i].im=out2[i][1]/(Num);
    }
    
}

void Compensator::write_init_and_difr()
{
    FILE* file1 = fopen("init_signal.txt","w");
    FILE* file2 = fopen("difr_signal.txt","w");
    for(int i =0 ; i<Num; i++)
    {
        fprintf(file1,"%f\n",P[i]);
        fprintf(file2,"%f\n",P3[i].re);
    }
      

}

void Compensator::init_signal_and_spectrum()
{
    for(int i=0; i<Num; i++)
    {
        t[i]=i*dt;//в секундах
        f_[i]=i*df;//в герцах
        P[i]=exp(-pow(t[i]-T/2,2)/pow(tau_0/1000000000,2));//tau_0 в наносекундах просто гаус
    }

    fftwf_plan plan;
    for (int i=0;i<Num;i++)
    {
        in[i][0]=P[i];
        in[i][1]=0;
    }
    plan = fftwf_plan_dft_1d(Num, in, out,FFTW_FORWARD,FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    for(int i=0;i<Num;i++)
    {
        P_sp[i].re=out[i][0];//init specrtum
        P_sp[i].im=out[i][1];
    }

}




/*
void Compensator::recovery_signal()
{
    //P3 - дифрагированный сигнал с учетом размера применика
    for (int i=0; i < Num; i++)
    {
        in[i][0]=P3[i].re;
        in[i][1] =0;
    }
    fftwf_plan plan1;
    plan1 = fftwf_plan_dft_1d(Num, in, out,FFTW_FORWARD,FFTW_ESTIMATE);
    fftwf_execute(plan1);
    fftwf_destroy_plan(plan1);
    for (int i=0; i < Num; i++)
    {
        P_sp[i].re = out[i][0];
        P_sp[i].im = out[i][1];
    }
    for (int i=1; i<Num;i++)//домножаем спектр нулевую частоту пропустил
    {

        P_rec[i]=recovery_one_f(P_sp[i],f_[i]);
        in[i][0]=P_rec[i].re;
        in[i][1]=P_rec[i].im;
    }
    fftwf_plan plan2;
    plan2 = fftwf_plan_dft_1d(Num, in, out,FFTW_BACKWARD,FFTW_ESTIMATE);
    fftwf_execute(plan2);
    fftwf_destroy_plan(plan2);
    for(int i=0;i<Num;i++)
    {
        P_rec[i].re=out[i][0]/(Num);
        P_rec[i].im=out[i][1]/(Num);
    }

}*/
/*void Compensator::draw_recovery_signal()
{
    double max = P_rec[0].re;
    int max_n = 0;
    for (int i =0; i<Num;i++)
    {
        if (P_rec[i].re >= max)
        {max = P_rec[i].re; max_n =i;}
    }

}*/

/*
void Compensator::on_init_signal_and_spectrum_clicked()
{
    Graph();
    Init_signal_and_spectrum();
    Graph3();
    F_d_const = c*1000*L/(pow(a,2)*M_PI);
}
*/



/*complex Compensator::recovery_one_f(complex P0,double f0)//на вход частота в мегагерцах
{

    if(rigth_compensation)
    {
        if (f0>F_N) {f0 = f0 - 2*F_N;}
        double B = pow(r0,2)/(a*a*(1+pow(cmp2(a,c,f0,L, f_d),2)));
        double ex1 = exp(-B);
        double module = M_PI * pow(a,2) * (pow(1.0 - ex1* cos(cmp2(a,c,f0,L, f_d)*B),2) + pow(ex1* sin(cmp2(a,c,f0,L, f_d)*B),2));
        complex q=complex((1.0 - ex1* cos(cmp2(a,c,f0,L, f_d)*B))/module,
                      -ex1* sin(cmp2(a,c,f0,L, f_d)*B)/module);
        return P0*q;
    }
    else
    {
        if (f0>F_N) {f0 = f0 - 2*F_N;}//когда этого нет а есть f-F_N в теле то работает

        complex q = complex(1.0  , -cmp2(a,c,f0,L, f_d));//получается что тут  1 - i*f_d/f
        return P0*q;

    }
}
*/
