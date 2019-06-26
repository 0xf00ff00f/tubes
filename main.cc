#include "panic.h"

#include "geometry.h"
#include "shader_program.h"
#include "noise.h"
#include "util.h"

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_inverse.hpp>

#include <algorithm>
#include <iostream>
#include <random>
#include <memory>

// #define DUMP_FRAMES

#define PERTURB

namespace
{
float rand_float()
{
    static std::mt19937 mt;
    return static_cast<float>(mt() - mt.min())/(mt.max() - mt.min());
}

int rand_int()
{
    static std::mt19937 mt;
    return mt();
}
}

class tube_geometry
{
public:
    tube_geometry(const glm::vec3 &p0, const glm::vec3 &p1)
        : p0_(p0)
        , p1_(0.5f*(p0 + p1))
        , p2_(p1)
    {
#ifdef PERTURB
        {
            glm::vec3 up(1, 1, 1);
            const auto dir = glm::normalize(p1 - p0);
            const auto side = glm::normalize(glm::cross(dir, up));
            up = glm::cross(side, dir);

            const float a = rand_float()*2.f*M_PI;
            const float x = sinf(a);
            const float y = cosf(a);
            p1_ += (side*x + up*y)*(0.5f + 1.5f*rand_float())*glm::length(p1 - p0);;
        }
#endif
        initialize_geometry();
    }

    void render() const
    {
        geometry_.bind();
        glDrawElements(GL_TRIANGLES, indices_.size(), GL_UNSIGNED_SHORT, nullptr);
    }

    glm::vec3 pos_at(float t) const
    {
        return p1_ + (1 - t)*(1 - t)*(p0_ - p1_) + t*t*(p2_ - p1_);
    }

    glm::vec3 dir_at(float t) const
    {
        return glm::normalize(2*(1 - t)*(p1_ - p0_) + 2*t*(p2_ - p1_));
    }

private:
    void initialize_geometry()
    {
        glm::vec3 up(1, 1, 1);

        for (int i = 0; i < rings_per_segment; ++i)
        {
            const auto t = static_cast<float>(i)/(rings_per_segment - 1);
            const auto center = pos_at(t);
            const auto dir = dir_at(t);
            const auto side = glm::normalize(glm::cross(dir, up));
            up = glm::cross(side, dir);

            const auto radius = ring_radius * (0.5f*(1.f - (0.75f*cosf(t*2.f*M_PI - M_PI) + 0.25f)) + 0.5f);

            for (int j = 0; j < verts_per_ring; ++j)
            {
                const auto a = static_cast<float>(j)*2.f*M_PI/verts_per_ring;
                const auto x = radius*cosf(a);
                const auto y = radius*sinf(a);
                const auto p = center + side*x + up*y;
                const auto n = glm::normalize(p - center);
                verts_.push_back({p, n});
            }
        }

        const auto total_rings = static_cast<int>(verts_.size())/verts_per_ring;
        for (int i = 0; i < total_rings - 1; ++i)
        {
            for (int j = 0; j < verts_per_ring; ++j)
            {
                const auto v00 = i*verts_per_ring + j;
                const auto v01 = i*verts_per_ring + (j + 1)%verts_per_ring;

                const auto v10 = (i + 1)*verts_per_ring + j;
                const auto v11 = (i + 1)*verts_per_ring + (j + 1)%verts_per_ring;

                indices_.push_back(v00);
                indices_.push_back(v01);
                indices_.push_back(v11);

                indices_.push_back(v11);
                indices_.push_back(v10);
                indices_.push_back(v00);
            }
        }

        geometry_.set_data(verts_, indices_);
    }

    static constexpr auto rings_per_segment = 40;
    static constexpr auto verts_per_ring = 30;
    static constexpr auto ring_radius = .15f;

    glm::vec3 p0_;
    glm::vec3 p1_;
    glm::vec3 p2_;
    using vertex = std::tuple<glm::vec3, glm::vec3>;
    std::vector<vertex> verts_;
    std::vector<GLshort> indices_;
    geometry geometry_;
};

constexpr auto window_width = 512;
constexpr auto window_height = 512;

int main(int argc, char *argv[])
{
    if (!glfwInit())
        panic("glfwInit failed\n");

    glfwSetErrorCallback([](int error, const char *description) { panic("GLFW error: %s\n", description); });

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_SAMPLES, 16);
    auto *window = glfwCreateWindow(window_width, window_height, "demo", nullptr, nullptr);
    if (!window)
        panic("glfwCreateWindow failed\n");

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    glewInit();

    glfwSetKeyCallback(window, [](GLFWwindow *window, int key, int scancode, int action, int mode) {
        if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
            glfwSetWindowShouldClose(window, GL_TRUE);
    });

    const auto projection =
        glm::perspective(glm::radians(45.0f), static_cast<float>(window_width) / window_height, 0.1f, 100.f);

    shader_program program;
    program.add_shader(GL_VERTEX_SHADER, "shaders/tube.vert");
    program.add_shader(GL_FRAGMENT_SHADER, "shaders/tube.frag");
    program.link();

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glEnable(GL_MULTISAMPLE);

#ifdef DUMP_FRAMES
    constexpr auto total_frames = 150;
    auto frame_num = 0;
#endif
    float cur_time = 0.0;

    {
    constexpr auto grid_cols = 5;
    constexpr auto grid_rows = 5;
    constexpr auto grid_layers = 5;

    constexpr auto cell_size = 5.;

    constexpr auto grid_width = cell_size * (grid_cols - 1);
    constexpr auto grid_height = cell_size * (grid_rows - 1);
    constexpr auto grid_depth = cell_size * (grid_layers - 1);

    std::vector<std::vector<std::vector<glm::vec3>>> grid(grid_layers);
    for (auto &layer : grid)
    {
        layer.resize(grid_rows);
        for (auto &row : layer)
            row.resize(grid_cols);
    }

    for (int i = 0; i < grid_layers; ++i)
    {
        for (int j = 0; j < grid_rows; ++j)
        {
            for (int k = 0; k < grid_cols; ++k)
            {
#ifdef PERTURB
                const auto z = -0.5*grid_depth + (i + (rand_float() - 0.5f))*cell_size;
                const auto y = -0.5*grid_height + (j + (rand_float() - 0.5f))*cell_size;
                const auto x = -0.5*grid_width + (k + (rand_float() - 0.5f))*cell_size;
#else
                const auto z = -0.5*grid_depth + i*cell_size;
                const auto y = -0.5*grid_height + j*cell_size;
                const auto x = -0.5*grid_width + k*cell_size;
#endif
                grid[i][j][k] = glm::vec3(x, y, z);
            }
        }
    }

    std::vector<std::unique_ptr<tube_geometry>> tubes;
    for (int i = 0; i < grid_layers; ++i)
    {
        for (int j = 0; j < grid_rows; ++j)
        {
            for (int k = 0; k < grid_cols; ++k)
            {
                if (i < grid_layers - 1)
                    tubes.emplace_back(new tube_geometry(grid[i][j][k], grid[i + 1][j][k]));

                if (j < grid_rows - 1)
                    tubes.emplace_back(new tube_geometry(grid[i][j][k], grid[i][j + 1][k]));

                if (k < grid_cols - 1)
                    tubes.emplace_back(new tube_geometry(grid[i][j][k], grid[i][j][k + 1]));
            }
        }
    }

    constexpr int num_lights = 17;
    struct light
    {
        tube_geometry *tube;
        float t;
        bool direction;
    };
    std::array<light, num_lights> lights;
    for (auto &light : lights)
    {
        light.tube = tubes[rand_int() % tubes.size()].get();
        light.t = rand_float();
        light.direction = rand_int() & 1;
    }

    while (!glfwWindowShouldClose(window))
    {
        glViewport(0, 0, window_width, window_height);
        glClearColor(0, 0, 0, 0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

#ifdef DUMP_FRAMES
        float angle = 1.5f * sinf(static_cast<float>(frame_num) * 2.0f * M_PI / total_frames);
#else
        float angle = 1.5f * sinf(cur_time);
#endif

        const auto view_pos = glm::vec3(0, 0, 20.f - 1.f*cur_time);
        const auto view_up = glm::vec3(cosf(0.1f*cur_time), sinf(0.1f*cur_time), 0);
        const auto view = glm::lookAt(view_pos, glm::vec3(0, 0, 0), view_up);

        // const auto model = glm::rotate(glm::mat4(1.0f), angle, glm::vec3(0, 1, 0));
        const auto model = glm::mat4(1.0f);
        const auto mvp = projection * view * model;

        glm::mat3 model_normal(model);
        model_normal = glm::inverse(model_normal);
        model_normal = glm::transpose(model_normal);

        program.bind();
        program.set_uniform(program.uniform_location("mvp"), mvp);
        program.set_uniform(program.uniform_location("normalMatrix"), model_normal);
        program.set_uniform(program.uniform_location("modelMatrix"), model);

        glm::vec3 light_position = tubes.front()->pos_at(glm::fract(cur_time));

        std::vector<glm::vec3> light_positions;
        std::transform(lights.begin(), lights.end(), std::back_inserter(light_positions), [](auto &light) {
            return light.tube->pos_at(light.t);
        });
        program.set_uniform(program.uniform_location("light_positions"), light_positions);
        program.set_uniform(program.uniform_location("global_light"), glm::vec3(5, 7, 5));

        for (auto &tube : tubes)
        {
            tube->render();
        }

        for (auto &light : lights)
        {
            constexpr float dt = 0.05;
            if (light.direction)
                light.t -= dt;
            else
                light.t += dt;
            if (light.t < 0 || light.t > 1)
            {
                light.tube = tubes[rand_int() % tubes.size()].get();
                light.direction = rand_int() & 1;
                light.t = light.direction ? 1 : 0;
            }
        }

#ifdef DUMP_FRAMES
        char path[80];
        std::sprintf(path, "%05d.ppm", frame_num);
        dump_frame_to_file(path, window_width, window_height);

        if (++frame_num == total_frames)
            break;
        cur_time += 1.0f / 40;
#else
        cur_time += 1.0f / 60;
#endif

        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    }

    glfwDestroyWindow(window);
    glfwTerminate();
}
