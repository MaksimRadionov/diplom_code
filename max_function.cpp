//я считал что входной вектор устроен так:
//например сэмплы в точке с координатами current_width и current_height находятся в интервале:
//vector<double> samples[((current_height - 1) * width + current_width - 1) * 4097 ....... {} * 4097 + 4096 ]
//current_height изменяется [1 .. height]
//в одной выборке 4097 точек, так как первый байт - синхронизация
//по обозначениям X -> current_width Y->current_height
//
//
#include <vector>
#include <algorithm>
#include <tuple>
double by_3_points(std::vector<double>& y, int n)
{
    return n + (y[n-1] - y[n+1])/(2 *(y[n-1] - 2  * y[n] + y[n+1]));
}

std::vector<std::tuple<double, int, int>> get_velocity(std::vector<double>& samples, float x1, float x2,// - я думаю тут вектор по ссылке идет
                float y1, float y2, int surf_z1, int surf_z2, int bottom_z1, int bottom_z2, 
                float height_z, bool max1, bool max2, float step_x, float step_y, int height, int width)//если расматривать ветор, как написано
                                                                                                     //то высота не нужна.       
{
    int down_width = int(x1 / step_x);
    int up_width = int(ceil(x2 / step_x));
    int down_height = int(y1 / step_y);
    int up_height = int(ceil(y2 / step_y));
    std::vector<std::tuple<double, int, int>> out;
    for(int current_height = down_height; current_height <= up_height; current_height++)
    {
        for(int current_width = down_width; current_width <= up_width; current_width++)
        { 
            std::vector<double>::iterator max_n; 
            if(max1)
                max_n = 
                        std::max_element(samples.begin() + ((current_height - 1) * width + current_width - 1) * 4097 + surf_z1,
                                samples.begin() + ((current_height - 1) * width + current_width - 1) * 4097 + surf_z2 );
            else
                max_n = 
                        std::min_element(samples.begin() + ((current_height - 1) * width + current_width - 1) * 4097 + surf_z1,
                                samples.begin() + ((current_height - 1) * width + current_width - 1) * 4097 + surf_z2 );
            int time_surf_n = max_n - samples.begin();
            double time_surf = by_3_points(samples, time_surf_n);    
        
            if(max2)
                max_n = 
                        std::max_element(samples.begin() + ((current_height - 1) * width + current_width - 1) * 4097 + bottom_z1,
                                samples.begin() + ((current_height - 1) * width + current_width - 1) * 4097 + bottom_z2 );
            else
                max_n = 
                        std::min_element(samples.begin() + ((current_height - 1) * width + current_width - 1) * 4097 + bottom_z1,
                                samples.begin() + ((current_height - 1) * width + current_width - 1) * 4097 + bottom_z2 );
            int time_bottom_n = max_n - samples.begin();
            double time_bottom = by_3_points(samples, time_bottom_n);

            out.push_back(std::make_tuple(2.0 * height_z / (time_bottom - time_surf),current_width, current_height));

        }
    }
    return out;
}  

