#include <sil/sil.hpp>
#include "random.hpp"
#include <cmath>
#include <complex>
#include <glm/gtx/matrix_transform_2d.hpp>
#include <vector>
#include <algorithm>

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

void change_brightness(sil::Image& image, float increase = 1.f){ // ⭐⭐ Luminosité
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

void convolution_box_blur(sil::Image& image, int offset = 1){ // ⭐⭐ Filtres séparables
    convolution(image, create_row_blur_kernel(offset));
    convolution(image, create_column_blur_kernel(offset));
}

void box_difference(sil::Image& image, float T = 25.f, float treshold = 0.1f, int blur1 = 3, int blur2 = 4){ // ⭐⭐ Différence de gaussiennes (de box_blur en l'occurence)
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

bool elements_in(std::vector<glm::vec3> const& v1, std::vector<glm::vec3> const& v2){ // Vérifie si tous les éléments de v1 sont dans v2
    int n = v1.size();
    int m = v2.size();
    if (n > m)
        return false;

    for (int i{0}; i < n; i++){
        bool is_in{false};
        glm::vec3 elem{v1[i]};
        for (int j{0}; j < m; j++){
            if (elem == v2[j]){
                is_in = true;
                break;
            }
        }
        if (!is_in)
            return false;
    }
    return true;
}

bool same_elements(std::vector<glm::vec3> const& v1, std::vector<glm::vec3> const& v2){ // Vérifie si v1 et v2 ont les mêmes éléments
    if (v1.size() != v2.size())
        return false;

    if (!elements_in(v1, v2))
        return false;

    return elements_in(v2, v1);
}

std::vector<glm::vec3> k_means(std::vector<glm::vec3> const& v, int k = 2, int max_iter = 100){ // Trouve les k couleurs les plus représentatives dans v
    std::vector<glm::vec3> means(k, glm::vec3{0.f});
    std::vector<glm::vec3> means_new(k, glm::vec3{0.f});
    std::vector<glm::vec3> sums(k, glm::vec3{0.f});
    std::vector<float> numbers(k, 0.f);

    for (int i{0}; i < k; i++){ // Initialisation des moyennes (aléatoire)
        means[i] = glm::vec3{random_float(0.f, 1.f)};
    }

    for (int i{0}; i < max_iter; i++){
        means_new = std::vector<glm::vec3>(k, glm::vec3{0.f});
        numbers = std::vector<float>(k, 0.f);

        // On associe chaque element de v à une moyenne et on les ajoute aux sommes
        for (int j{0}; j < v.size(); j++){
            glm::vec3 color_j{v[j]};

            // On cherche la moyenne la plus proche
            int nearest_l = 0;
            float nearest_distance = glm::distance(color_j, means[0]);
            for (int l{1}; l < k; l++){
                float distance = glm::distance(color_j, means[l]);
                if (distance < nearest_distance){
                    nearest_l = l;
                    nearest_distance = distance;
                }
            }

            //On ajoute la couleur
            means_new[nearest_l] += color_j;
            numbers[nearest_l]++;
        }

        // On calcule chaque nouvelle moyenne
        for (int i{0}; i < k; i++){
            if (numbers[i] > 0)
                means_new[i] /= numbers[i];
            else
                means_new[i] = glm::vec3{random_float(0.f, 1.f)};
        }

        if (same_elements(means, means_new)){
            //std::cout << "Iterations : " << i + 1 << std::endl;
            break;
        }

        means = means_new;
    }

    return means;
}

void set_colors(sil::Image& image, std::vector<glm::vec3> v){ // Transforme les couleurs de l'image en la couleur la plus proche dans v
    for (int x{0}; x < image.width(); x++)
    {
        for (int y{0}; y < image.height(); y++)
        {
            glm::vec3 color = image.pixel(x, y);

            // On cherche la couleur
            int nearest_l = 0;
            float nearest_distance = glm::distance(color, v[0]);
            for (int l{1}; l < v.size(); l++){
                float distance = glm::distance(color, v[l]);
                if (distance < nearest_distance){
                    nearest_l = l;
                    nearest_distance = distance;
                }
            }

            image.pixel(x, y) = v[nearest_l];
        }
    }
}

void k_colors(sil::Image& image, int k = 2, int max_iter = 100){ // ⭐⭐⭐⭐⭐ K-means : trouver les couleurs les plus présentes dans une image
    set_colors(image, k_means(image.pixels(), k, max_iter));
}

template <typename T> int index_min(std::vector<T> v){
    int i_min{0};
    T val_min{v[0]};
    for (int i{1}; i < v.size(); i++){
        T val = v[i];
        if (val < val_min){
            i_min = i;
            val_min = val;
        }
    }
    return i_min;
}

glm::vec3 vec3_mean(std::vector<glm::vec3> v){
    glm::vec3 v_sum{0.f};
    for (glm::vec3 val : v){
        v_sum += val;
    }
    return v_sum/static_cast<float>(v.size());
}

float vec3_std(std::vector<glm::vec3> v, glm::vec3 v_mean){
    float v_std{0.f};
    for (glm::vec3 val : v){
        v_std += glm::distance(val, v_mean);
    }
    return v_std/static_cast<float>(v.size());
}

void statistics_square(sil::Image const& image, int left_x, int up_y, int side_size, std::vector<glm::vec3>& means, std::vector<float>& st_deviations){
    int limit_x{image.width()}, limit_y{image.height()};
    std::vector<glm::vec3> v{};

    for (int x{left_x}; x < left_x + side_size; x++)
        for (int y{up_y}; y < up_y + side_size; y++){
            int x2{x}, y2{y};

            if (x2 < 0)
                x2 = 0;
            else if (limit_x <= x2)
                x2 = limit_x - 1;

            if (y2 < 0)
                y2 = 0;
            else if (limit_y <= y2)
                y2 = limit_y - 1;

            v.push_back(image.pixel(x2, y2));
        }

    glm::vec3 v_mean = vec3_mean(v);
    means.push_back(v_mean);

    st_deviations.push_back(vec3_std(v, v_mean));
}

glm::vec3 kuwahara_pixel(sil::Image const& image, int offset, int x, int y){
    std::vector<glm::vec3> means{};
    std::vector<float> st_deviations{};
    
    statistics_square(image, x - offset, y - offset, offset + 1, means, st_deviations); // Haut gauche
    statistics_square(image, x,          y - offset, offset + 1, means, st_deviations); // Haut droite
    statistics_square(image, x - offset, y,          offset + 1, means, st_deviations); // Bas  gauche
    statistics_square(image, x,          y,          offset + 1, means, st_deviations); // Bas  droite
    
    return means[index_min(st_deviations)];
}

void kuwahara_filter(sil::Image& image, int offset = 2){ // ⭐⭐⭐⭐⭐ Filtre de Kuwahara (effet peinture à l'huile) 28
    sil::Image image_copy{image};
    for (int x{0}; x < image.width(); x++)
    {
        for (int y{0}; y < image.height(); y++)
        {
            image.pixel(x, y) = kuwahara_pixel(image_copy, offset, x, y);
        }
    }
}

void diamond_step(sil::Image& image, int half_size, int x, int y, float rand_range){
    float sum{0.f};
    int n = 0;

    bool affiche{x == 128 && y == 256};

    if (0 <= x - half_size){
        sum += image.pixel(x - half_size, y).r;
        n++;
        if (affiche) std::cout << "pixel1=" << image.pixel(x - half_size, y).r << std::endl;
    }
    if (x + half_size < image.width()){
        sum += image.pixel(x + half_size, y).r;
        n++;
        if (affiche) std::cout << "pixel2=" << image.pixel(x + half_size, y).r << std::endl;
    }
    if (0 <= y - half_size){
        sum += image.pixel(x, y - half_size).r;
        n++;
        if (affiche) std::cout << "pixel3=" << image.pixel(x, y - half_size).r << std::endl;
    }
    if (y + half_size < image.height()){
        sum += image.pixel(x, y + half_size).r;
        n++;
        if (affiche) std::cout << "pixel4=" << image.pixel(x, y + half_size).r << std::endl;
    }

    if (affiche){
        std::cout << "sum=" << sum << ", n=" << n << ", res=" << sum/static_cast<float>(n) << std::endl;

    }

    image.pixel(x, y) = glm::vec3{sum/static_cast<float>(n) + random_float(-rand_range, rand_range)};
}

sil::Image diamond_square(int n = 9, float min_h = 0.f, float max_h = 1.f, float rand_range = 0.5f){ // ⭐⭐⭐⭐⭐⭐ Diamond Square 29
    int square_size = pow(2, n); // Taille d'un carré - 1
    int map_size = square_size + 1;
    int n_squares = 1; // Nombre de carrés par ligne

    sil::Image image = sil::Image{map_size, map_size};

    // Initialisation des coins
    image.pixel(0, 0) = glm::vec3{random_float(min_h, max_h)};
    image.pixel(square_size, 0) = glm::vec3{random_float(min_h, max_h)};
    image.pixel(0, square_size) = glm::vec3{random_float(min_h, max_h)};
    image.pixel(square_size, square_size) = glm::vec3{random_float(min_h, max_h)};

    while (square_size > 1){
        //int x_first = square_size/2; // Abcisse du centre des carrés de gauche
        int half_size = square_size/2;

        for (int i = 0; i < n_squares; i++){
            int x_left = i * square_size; 
            int x_right = x_left + square_size;
            int x_center = x_left + half_size;

            for (int j = 0; j < n_squares; j++){
                int y_up = j * square_size; 
                int y_down = y_up + square_size;
                int y_center = y_up + half_size;
                
                image.pixel(x_center, y_center) = glm::vec3{random_float(-rand_range, rand_range)} + (
                    image.pixel(x_left, y_up) +
                    image.pixel(x_right, y_up) + 
                    image.pixel(x_left, y_down) + 
                    image.pixel(x_right, y_down))/4.f;
            }
        }

        for (int i = 0; i < n_squares; i++){
            int x_left = i * square_size; 
            int x_right = x_left + square_size;
            int x_center = x_left + half_size;

            for (int j = 0; j < n_squares; j++){
                int y_up = j * square_size; 
                int y_down = y_up + square_size;
                int y_center = y_up + half_size;
                
                if (i == 0)
                    diamond_step(image, half_size, x_left, y_center, rand_range); // Gauche
                if (j == 0)
                    diamond_step(image, half_size, x_center, y_up, rand_range); // Haut
                diamond_step(image, half_size, x_right, y_center, rand_range); // Droite
                diamond_step(image, half_size, x_center, y_down, rand_range); // Bas
            }
        }

        rand_range /= 2;
        n_squares *= 2;
        square_size = square_size/2;
    }

    return image;
}

float round_to(float x, int decimals){
    int coeff = pow(10, decimals);
    return std::round(x * coeff)/coeff;
}

void color_height_map(sil::Image& image, float treshold = 0.5f){ // ⭐⭐ Colorer la height map
    for (glm::vec3& color : image.pixels())
    {
        float pixel_height{color.r};
        color.r = 0.f;
        if (pixel_height > treshold){
            color.g = 1.f - pixel_height;
            color.b = 0.f;
        }
        else{
            color.g = 0.f;
        }
        color.r = round_to(color.r, 1);
        color.g = round_to(color.g, 1);
        color.b = round_to(color.b, 1);
    }
}

int main()
{
    std::cout << "Preparation" << std::endl;

    //std::vector<std::vector<float>> kernel_box_blur_10{create_box_blur_kernel(10)};
    //std::vector<std::vector<float>> kernel_box_blur_50{create_box_blur_kernel(50)};
    
    //sil::Image image{"images/logo.png"};
    //sil::Image image{"images/photo_faible_contraste.jpg"};
    //sil::Image image{"images/photo.jpg"};
    sil::Image image{"output/diamond_square.png"};

    std::cout << "Execution" << std::endl;

    color_height_map(image);
    //image = diamond_square();

    image.save("output/color_height_map.png");
    image.save("output/pouet.png");

    std::cout << "Fin" << std::endl;
}
