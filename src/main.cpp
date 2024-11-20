#include <sil/sil.hpp>
#include "random.hpp"
#include <cmath>
#include <complex>
#include <glm/gtx/matrix_transform_2d.hpp>

#include <vector>
#include <string>
#include <iostream>

#define PI 3.14159265

void keep_green_only(sil::Image& image){ // ⭐ Ne garder que le vert
    for (glm::vec3& color : image.pixels())
    {
        color.r = color.b = 0.f;
    }
}

void channels_swap(sil::Image& image){ // ⭐ Échanger les canaux
    float temp = 0.f;
    for (glm::vec3& color : image.pixels())
    {
        temp = color.r;
        color.r = color.b;
        color.b = temp;
    }
}

void black_and_white(sil::Image& image){ // ⭐ Noir & blanc
    for (glm::vec3& color : image.pixels())
    {
        float luminance = 0.2126f * color.r + 0.7152f * color.g + 0.0722f * color.b;
        color.r = color.g = color.b = luminance;
    }
}

void negative(sil::Image& image){ // ⭐ Négatif
    for (glm::vec3& color : image.pixels())
    {
        color = 1.f - color;
    }
}

sil::Image gradient(int width = 300, int height = 200){ // ⭐ Dégradé
    sil::Image image{width, height};
    for (int x{0}; x < width; x++)
    {
        for (int y{0}; y < height; y++)
        {
            image.pixel(x, y) = glm::vec3{static_cast<float>(x)/(static_cast<float>(width) - 1)};
        }
    }
    return image;
}

void miroir(sil::Image& image){ // ⭐⭐ Miroir
    sil::Image image_copy{image};
    for (int x{0}; x < image.width(); x++)
    {
        for (int y{0}; y < image.height(); y++)
        {
            image.pixel(x, y) = image_copy.pixel(image.width() - x - 1, y);
        }
    }
}

void noise(sil::Image& image, float noise_level = 0.3f){ // ⭐⭐ Image bruitée
    sil::Image image_copy{image};
    for (glm::vec3& color : image.pixels())
    {
        if (true_with_probability(noise_level)){
            color = {random_float(0.f, 1.f), random_float(0.f, 1.f), random_float(0.f, 1.f)};
        }
    }
}

void rotation_90(sil::Image& image){ // ⭐⭐ Rotation de 90°
    sil::Image image_res{image.height(), image.width()};
    for (int x{0}; x < image_res.width(); x++)
    {
        for (int y{0}; y < image_res.height(); y++)
        {
            image_res.pixel(x, y) = image.pixel(y, image.height() - x -1);
        }
    }
    image = image_res;
}

void rgb_split(sil::Image& image, int decalage = 30){ // ⭐⭐ RGB split
    sil::Image image_copy{image};
    for (int x{0}; x < image.width(); x++)
    {
        for (int y{0}; y < image.height(); y++)
        {
            if (x + decalage < image.width())
                image.pixel(x, y).b = image_copy.pixel(x + decalage, y).g;
            else
                image.pixel(x, y).b = 0.f;
            
            if (x - decalage >= 0)
                image.pixel(x, y).r = image_copy.pixel(x - decalage, y).r;
            else
                image.pixel(x, y).r = 0.f;
        }
    }
}

void change_brightness(sil::Image& image, int increase = 1){ // ⭐⭐ Luminosité
    // Plus increase est grand, plus l'image est lumineuse
    // Plus increase est loin en dessous de 0, plus l'image est sombre
    float power = pow(2, -increase);
    for (glm::vec3& color : image.pixels())
    {
        color.r = pow(color.r, power);
        color.g = pow(color.g, power);
        color.b = pow(color.b, power);
    }
}

sil::Image disc(float radius = 100.f, float center_x = 249.5f, float center_y = 249.5f, int width = 500, int height = 500){ // ⭐⭐(⭐) Disque
    sil::Image image{width, height};
    float distance = 0.f;
    for (int x{0}; x < width; x++)
    {
        for (int y{0}; y < height; y++)
        {
            distance = pow(pow(x - center_x, 2) + pow(y - center_y, 2), 0.5f);
            if (distance <= radius)
                image.pixel(x, y) = glm::vec3{1.f};
        }
    }
    return image;
}

void add_circle(sil::Image& image, float radius, float thickness, float center_x, float center_y, int width, int height){ // Pour utiliser dans les fonctions suivantes
    float distance = 0.f;
    for (int x{0}; x < width; x++)
    {
        for (int y{0}; y < height; y++)
        {
            distance = pow(pow(x - center_x, 2) + pow(y - center_y, 2), 0.5f);
            if (std::abs(distance - radius) <= thickness/2.f)
                image.pixel(x, y) = glm::vec3{1.f};
        }
    }
}

sil::Image circle(float radius = 100.f, float thickness = 10.f, float center_x = 249.5f, float center_y = 249.5f, int width = 500, int height = 500){ // ⭐ Cercle
    sil::Image image{width, height};
    add_circle(image, radius, thickness, center_x, center_y, width, height);
    return image;
}

void animation(){ // ⭐⭐ Animation
    sil::Image image{"images/photo.jpg"};
    int n_frames = 40;
    float ecart = 500.f/(n_frames - 1);
    for(int i = 0; i<n_frames; i++){
        image = disc(100.f, 100.f, ecart * i);
        image.save("output/animation/cercle" + std::to_string(i) + ".png");
    }
}

sil::Image rosace(int n_petals = 6, float radius = 100.f, float thickness = 5.f, float center_x = 249.5f, float center_y = 249.5f, int width = 500, int height = 500){ // ⭐⭐⭐ Rosace
    sil::Image image{width, height};

    add_circle(image, radius, thickness, center_x, center_y, width, height);

    float offset_x{0.f}, offset_y{0.f};
    float angle{0.f};
    for(int i{0}; i<n_petals; i++){
        angle = (360.f/n_petals) * i;
        offset_x = std::cos(angle * PI / 180.0) * radius;
        offset_y = std::sin(angle * PI / 180.0) * radius;
        add_circle(image, radius, thickness, center_x + offset_x, center_y + offset_y, width, height);
    }

    return image;
}

void mosaic(sil::Image& image, int n_divisions = 5){ // ⭐⭐ Mosaïque
    sil::Image image_res{image.width(), image.height()};
    int division_width = image.width()/n_divisions;
    int division_height = image.height()/n_divisions;
    int p_x{0}, p_y{0};

    for (int x{0}; x < image_res.width(); x++)
    {
        for (int y{0}; y < image_res.height(); y++)
        {
            p_x = (x%division_width)*n_divisions;
            p_y = (y%division_height)*n_divisions;
            image_res.pixel(x, y) = image.pixel(p_x, p_y);
        }
    }

    image = image_res;
}

void mosaic_mirror(sil::Image& image, int n_divisions = 5){ // ⭐⭐⭐⭐ Mosaïque miroir
    sil::Image image_res{image.width(), image.height()};
    int division_width = image.width()/n_divisions;
    int division_height = image.height()/n_divisions;
    int p_x{0}, p_y{0};

    for (int x{0}; x < image_res.width(); x++)
    {
        for (int y{0}; y < image_res.height(); y++)
        {
            if ((x/division_width)%2 == 0)
                p_x = (x%division_width)*n_divisions;
            else
                p_x = ((image.width() - x - 1)%division_width)*n_divisions;
            
            if ((y/division_height)%2 == 0)
                p_y = (y%division_height)*n_divisions;
            else
                p_y = ((image.height() - y - 1)%division_height)*n_divisions;
            
            image_res.pixel(x, y) = image.pixel(p_x, p_y);
        }
    }

    image = image_res;
}

void swap_rectangles(sil::Image& image, int rect_width, int rect_height, int x1, int y1, int x2, int y2){ // Fonction pour utiliser dans glitch
    sil::Image image_copy{image};
    for (int dx{0}; dx < rect_width; dx++)
    {
        for (int dy{0}; dy < rect_height ; dy++)
        {
            image.pixel(x1 + dx, y1 + dy) = image_copy.pixel(x2 + dx, y2 + dy);
            image.pixel(x2 + dx, y2 + dy) = image_copy.pixel(x1 + dx, y1 + dy);
        }
    }
}

void glitch(sil::Image& image, int n_glitches = 100, int max_width = 30, int max_height = 10){ // ⭐⭐⭐ Glitch
    if (max_width>image.width())
        max_width = image.width();
    if (max_height>image.height())
        max_height = image.height();

    for(int i{0}; i<n_glitches; i++){
        swap_rectangles(image, random_int(0, max_width), random_int(0, max_height), random_int(0, image.width() - max_width), random_int(0, image.height() - max_height), random_int(0, image.width() - max_width), random_int(0, image.height() - max_height));
    }
}

float brightness(glm::vec3 const& color){ // Luminance d'une couleur
    return 0.2126f * color.r + 0.7152f * color.g + 0.0722f * color.b;
}

void sort_pixels(sil::Image& image){ // ⭐⭐⭐ Tri de pixels 18
    std::vector<glm::vec3> v = image.pixels();
    std::sort(image.pixels().begin(), image.pixels().end(), [](glm::vec3 const& color1, glm::vec3 const& color2)
    {
        return brightness(color1) < brightness(color2);
    });
}

void sort_pixels_by_row(sil::Image& image){ // Tri de pixels par ligne
    std::vector<glm::vec3> v = image.pixels();
    for(int i{0}; i<image.height(); i++){
        std::sort(image.pixels().begin() + i * image.width(), image.pixels().begin() + (i + 1) * image.width(), [](glm::vec3 const& color1, glm::vec3 const& color2)
        {
            return brightness(color1) < brightness(color2);
        });
    }
}

sil::Image gradient_colors(glm::vec3 color1 = {0.f, 1.f, 0.f}, glm::vec3 color2 = {1.f, 0.f, 0.f}, int width = 300, int height = 200){ // Dégradé de couleurs
    sil::Image image{width, height};
    float progression{0.f};
    for (int x{0}; x < width; x++)
    {
        for (int y{0}; y < height; y++)
        {
            progression = static_cast<float>(x)/(static_cast<float>(width) - 1);
            image.pixel(x, y) = glm::vec3{progression * color1 + (1.f - progression) * color2};
        }
    }
    return image;
}

// OKLAB A FAIRE

int iter_mandelbrot(std::complex<float> c, int iter_max = 10){
    std::complex<float> z{0.f, 0.f};
    int n = 0;
    while (std::abs(z) <= 2){
        n++;
        if (n>=iter_max)
            break;
        z = z * z + c;
    }
    return n;
}

sil::Image mandelbrot(int iter_max = 50, int width = 500, int height = 500){ // ⭐⭐⭐(⭐) Fractale de Mandelbrot
    sil::Image image{width, height};

    for (int x{0}; x < width; x++)
    {
        for (int y{0}; y < height; y++)
        {
            std::complex<float> z{x/(width - 1.f) * 4.f - 2.f, y/(height - 1.f) * 4.f - 2.f};
            image.pixel(x, y) = glm::vec3{iter_mandelbrot(z, iter_max)/static_cast<float>(iter_max)};
        }
    }
    
    return image;
}

float bayer_matrix_4x4[][4] = {
    {    -0.5,       0,  -0.375,   0.125 },
    {    0.25,   -0.25,   0.375, - 0.125 },
    { -0.3125,  0.1875, -0.4375,  0.0625 },
    {  0.4375, -0.0625,  0.3125, -0.1875 },
};

void dithering_4x4(sil::Image& image, float r = 1.f){ // ⭐⭐⭐(⭐) Tramage
    float gray_scale = 0.f;
    for (int x{0}; x < image.width(); x++)
    {
        for (int y{0}; y < image.height(); y++)
        {
            gray_scale = brightness(image.pixel(x, y)) + r * bayer_matrix_4x4[y%4][x%4];
            if (gray_scale > 0.5f)
                gray_scale = 1.f;
            else
                gray_scale = 0.f;
            image.pixel(x, y) = glm::vec3{gray_scale};
        }
    }
}

void new_min(float cur, float& min){
    if (cur < min)
        min = cur;
}

void new_max(float cur, float& max){
    if (cur > max)
        max = cur;
}

void bad_normalize(sil::Image& image){ // Normalisation de l'histogramme mais ça fait des couleurs bizarres
    glm::vec3 min_rgb{1.f, 1.f, 1.f};
    glm::vec3 max_rgb{0.f, 0.f, 0.f};
    glm::vec3 cur_rgb{0.f, 0.f, 0.f};

    for (int x{0}; x < image.width(); x++)
    {
        for (int y{0}; y < image.height(); y++)
        {
            cur_rgb = image.pixel(x, y);

            new_min(cur_rgb.r, min_rgb.r);
            new_min(cur_rgb.g, min_rgb.g);
            new_min(cur_rgb.b, min_rgb.b);

            new_max(cur_rgb.r, max_rgb.r);
            new_max(cur_rgb.g, max_rgb.g);
            new_max(cur_rgb.b, max_rgb.b);
        }
    }

    glm::vec3 dif_rgb = max_rgb - min_rgb;

    for (int x{0}; x < image.width(); x++)
    {
        for (int y{0}; y < image.height(); y++)
        {
            image.pixel(x, y) = (image.pixel(x, y) - min_rgb)/dif_rgb;
        }
    }
}

void normalize(sil::Image& image){ // ⭐⭐⭐(⭐) Normalisation de l'histogramme
    float min_bri{1.f};
    float max_bri{0.f};
    float cur_bri{0.f};

    for (int x{0}; x < image.width(); x++)
    {
        for (int y{0}; y < image.height(); y++)
        {
            cur_bri = brightness(image.pixel(x, y));
            new_min(cur_bri, min_bri);
            new_max(cur_bri, max_bri);
        }
    }

    float dif_bri = max_bri - min_bri;

    for (int x{0}; x < image.width(); x++)
    {
        for (int y{0}; y < image.height(); y++)
        {
            image.pixel(x, y) = (image.pixel(x, y) - min_bri)/dif_bri;
        }
    }
}

void modif(sil::Image& image){ // Fonction qui fait des couleurs bizarres (je voulais tester)
    glm::vec3 p{0.f};
    for (int x{0}; x < image.width(); x++)
    {
        for (int y{0}; y < image.height(); y++)
        {
            p = image.pixel(x, y);
            image.pixel(x, y).r = (0.f  ) * p.r + (0.f  ) * p.g + (1.f  ) * p.b;
            image.pixel(x, y).g = (1.f  ) * p.r + (0.f  ) * p.g + (0.f  ) * p.b;
            image.pixel(x, y).b = (0.f  ) * p.r + (1.f  ) * p.g + (0.f  ) * p.b;
        }
    }
}

glm::vec2 rotated(glm::vec2 point, glm::vec2 center_of_rotation, float angle)
{
    return glm::vec2{glm::rotate(glm::mat3{1.f}, angle) * glm::vec3{point - center_of_rotation, 0.f}} + center_of_rotation;
}

void vortex(sil::Image& image, float intensity = 0.1f){ // ⭐⭐⭐⭐ Vortex
    sil::Image image_copy{image};
    image = sil::Image{image.width(), image.height()};

    glm::vec2 coordinates{0.f, 0.f};
    glm::vec2 center{image.width()/2.f, image.height()/2.f};
    float distance = 0.f;
    float coor_x{0.f}, coor_y{0.f};

    for (int x{0}; x < image.width(); x++)
    {
        for (int y{0}; y < image.height(); y++)
        {
            coordinates = rotated(glm::vec2{x, y}, center, distance*intensity);
            distance = glm::distance(center, coordinates);
            coor_x = static_cast<int>(coordinates.x);
            coor_y = static_cast<int>(coordinates.y);
            if (0 <= coor_x && coor_x < image.width() && 0 <= coor_y && coor_y < image.height())
                image.pixel(x, y) = image_copy.pixel(coor_x, coor_y);
        }
    }
}

void convolution(sil::Image& image, std::vector<std::vector<float>> kernel){ // ⭐⭐⭐⭐ Convolutions
    sil::Image image_copy{image};
    glm::vec3 result{0.f};
    unsigned int n{kernel[0].size()}, m{kernel.size()};
    unsigned int offset_u{n/2}, offset_v{n/2};
    int x2{0}, y2{0};
    int width{image.width()}, height{image.height()};
    for (int x{0}; x < width; x++)
    {
        for (int y{0}; y < height; y++)
        {
            result = glm::vec3{0.f};
            for (int v{0}; v < m; v++){
                for (int u{0}; u < n; u++){
                    x2 = x + u - offset_u;
                    y2 = y + v - offset_v;

                    if (x2 < 0) 
                        x2 = 0;
                    else if (width <= x2) 
                        x2 = width - 1;

                    if (y2 < 0) 
                        y2 = 0;
                    else if (height <= y2) 
                        y2 = height - 1;

                    result += image_copy.pixel(x2, y2) * kernel[v][u];
                }
            }
            image.pixel(x, y) = result;
        }
    }
}

std::vector<std::vector<float>> create_rectangle_blur_kernel(int offset_u = 1, int offset_v = 1){
    int n = 1 + offset_u * 2;
    int m = 1 + offset_v * 2;
    std::vector<float> row(n, 1.f/(n * m));
    return std::vector<std::vector<float>>(m, row);
}

std::vector<std::vector<float>> create_box_blur_kernel(int offset = 1){
    return create_rectangle_blur_kernel(offset, offset);
}

// ⭐ Netteté, Contours, etc.

std::vector<std::vector<float>> kernel_emboss {{
    {-2.f, -1.f, 0.f},
    {-1.f, 1.f, 1.f},
    {0.f, 1.f, 2.f}}};

std::vector<std::vector<float>> kernel_outline {{
    {-1.f, -1.f, -1.f},
    {-1.f, 8.f, -1.f},
    {-1.f, -1.f, -1.f}}};

std::vector<std::vector<float>> kernel_sharpen {{
    {0.f, -1.f, 0.f},
    {-1.f, 5.f, -1.f},
    {0.f, -1.f, 0.f}}};

std::vector<std::vector<float>> create_row_blur_kernel(int offset = 1){
    return create_rectangle_blur_kernel(offset, 1);
}

std::vector<std::vector<float>> create_column_blur_kernel(int offset = 1){
    return create_rectangle_blur_kernel(1, offset);
}

void convolution_box_blur(sil::Image& image, int offset = 1){ // ⭐⭐ Filtres séparables 25
    convolution(image, create_row_blur_kernel(offset));
    convolution(image, create_column_blur_kernel(offset));
}

// (1 + T) * G1 - T * G2

void box_difference(sil::Image& image, float T = 25.f, float treshold = 0.1f, int blur1 = 3, int blur2 = 4){ // ⭐⭐ Différence de gaussiennes (de box_blur en l'occurence) 26
    sil::Image image_blur1{image};
    sil::Image image_blur2{image};

    convolution_box_blur(image_blur1, blur1);
    convolution_box_blur(image_blur2, blur2);

    float luminance1{0.f}, luminance2{0.f}, luminance_new{0.f};

    for (int x{0}; x < image.width(); x++)
    {
        for (int y{0}; y < image.height(); y++)
        {
            luminance1 = brightness(image_blur1.pixel(x, y));
            luminance2 = brightness(image_blur2.pixel(x, y));
            luminance_new = (T + 1) * luminance1 - T * luminance2;

            if (luminance_new > treshold)
                luminance_new = 1.f;
            else
                luminance_new = 0.f;

            image.pixel(x, y) = glm::vec3(luminance_new);
        }
    }

    normalize(image);
}

bool same_elements(std::vector<glm::vec3> v1, std::vector<glm::vec3> v2)
{
    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());
    return v1 == v2;
}

std::vector<glm::vec3> k_means(std::vector<glm::vec3> v, int k, int max_iter = 100){
    std::vector<glm::vec3> means{};
    std::vector<glm::vec3> means_new{};

    for (int i{0}; i < k; i++){
        means.push_back(v[i]);
    }

    for (int i{0}; i < max_iter; i++){
        for (int j{0}; j < max_iter; j++){
            
        }
    }
}

int main()
{
    std::cout<<"Preparation"<<std::endl;

    std::vector<std::vector<float>> kernel_box_blur_10{create_box_blur_kernel(10)};
    //std::vector<std::vector<float>> kernel_box_blur_50{create_box_blur_kernel(50)};
    
    //sil::Image image{"images/logo.png"};
    //sil::Image image{"images/photo_faible_contraste.jpg"};
    sil::Image image{"images/photo.jpg"};

    std::cout<<"Execution"<<std::endl;

    box_difference(image);
    //image = mandelbrot();

    image.save("output/pouet.png");

    std::cout<<"Fin"<<std::endl;
}

// 26/31
// OKLAB + gaussian blur, dithering général