#include "panic.h"

#include "geometry.h"
#include "shader_program.h"
#include "util.h"
#include "rand.h"
#include "blur_effect.h"

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_inverse.hpp>

#include <algorithm>
#include <iostream>
#include <memory>

// #define DUMP_FRAMES

#define PERTURB
#define GLOW

class sphere_geometry
{
public:
    sphere_geometry(glm::vec3 &center, float radius, int max_subdivisions)
        : center_(center)
        , radius_(radius)
        , max_subdivisions_(max_subdivisions)
    {
        initialize_geometry(center, radius);
    }

    void render() const
    {
        geometry_.bind();
        glDrawElements(GL_TRIANGLES, indices_.size(), GL_UNSIGNED_SHORT, nullptr);
    }

private:
    void initialize_geometry(const glm::vec3 &center, float radius)
    {
        const auto i1 = add_vertex(glm::vec3{0.000000, 0.000000, 0.612372});
        const auto i2 = add_vertex(glm::vec3{-0.288675, -0.500000, -0.204124});
        const auto i3 = add_vertex(glm::vec3{-0.288675, 0.500000, -0.204124});
        const auto i4 = add_vertex(glm::vec3{0.577350, 0.000000, -0.204124});

        subdivide_triangle(i2, i3, i4, 0);
        subdivide_triangle(i3, i2, i1, 0);
        subdivide_triangle(i4, i1, i2, 0);
        subdivide_triangle(i1, i4, i3, 0);

        geometry_.set_data(verts_, indices_);
    }

    void subdivide_triangle(int i0, int i1, int i2, int level)
    {
        if (level == max_subdivisions_)
        {
            indices_.push_back(i0);
            indices_.push_back(i1);
            indices_.push_back(i2);
        }
        else
        {
            const auto &v0 = std::get<1>(verts_[i0]);
            const auto &v1 = std::get<1>(verts_[i1]);
            const auto &v2 = std::get<1>(verts_[i2]);

            const auto i01 = add_vertex(0.5f*(v0 + v1));
            const auto i12 = add_vertex(0.5f*(v1 + v2));
            const auto i20 = add_vertex(0.5f*(v2 + v0));

            subdivide_triangle(i0, i01, i20, level + 1);
            subdivide_triangle(i01, i1, i12, level + 1);
            subdivide_triangle(i20, i12, i2, level + 1);
            subdivide_triangle(i01, i12, i20, level + 1);
        }
    }

    int add_vertex(const glm::vec3 &v)
    {
        const auto n = glm::normalize(v);
        const auto p = center_ + radius_*n;
        verts_.emplace_back(p, n);
        return verts_.size() - 1;
    }

    glm::vec3 center_;
    float radius_;
    int max_subdivisions_;

    using vertex = std::tuple<glm::vec3, glm::vec3>;
    std::vector<vertex> verts_;
    std::vector<GLshort> indices_;
    geometry geometry_;
};

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

        for (int i = 0; i < rings_per_segment - 1; ++i)
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

class demo
{
public:
    demo(int window_width, int window_height)
        : window_width_(window_width)
        , window_height_(window_height)
        , blur_(window_width/4, window_height/4)
    {
        initialize_shader();
        initialize_geometry();
    }

    void render_and_step(float dt)
    {
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LESS);

#ifdef GLOW
        blur_.bind();
        glViewport(0, 0, blur_.width(), blur_.height());
        glClearColor(0, 0, 0, 0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glDisable(GL_MULTISAMPLE);
        render(program_glow_);
#endif

        framebuffer::unbind();
        glViewport(0, 0, window_width_, window_height_);
        glClearColor(0, 0, 0, 0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glEnable(GL_MULTISAMPLE);
        render(program_);

#ifdef GLOW
        blur_.render(window_width_, window_height_);
#endif

        cur_time_ += dt;
        for (auto &light : lights_)
        {
            constexpr auto light_speed = 3.f;
            light.t += light.direction ? -dt*light_speed : dt*light_speed;
            if (light.t < 0 || light.t > 1)
            {
                light.tube = tubes_[rand_int() % tubes_.size()].get();
                light.direction = rand_int() & 1;
                light.t = light.direction ? 1 : 0;
            }
        }
    }

private:
    void initialize_shader()
    {
        program_.add_shader(GL_VERTEX_SHADER, "shaders/tube.vert");
        program_.add_shader(GL_FRAGMENT_SHADER, "shaders/tube.frag");
        program_.link();

        program_.bind();
        program_.set_uniform(program_.uniform_location("global_light"), glm::vec3(5, 7, 5));

        program_glow_.add_shader(GL_VERTEX_SHADER, "shaders/tube.vert");
        program_glow_.add_shader(GL_FRAGMENT_SHADER, "shaders/tube-glow.frag");
        program_glow_.link();
    }

    void initialize_geometry()
    {
        // grid

        for (int i = 0; i < grid_layers; ++i)
        {
            for (int j = 0; j < grid_rows; ++j)
            {
                for (int k = 0; k < grid_cols; ++k)
                {
                    const auto z = -0.5f*grid_depth + i*cell_size;
                    const auto y = -0.5f*grid_height + j*cell_size;
                    const auto x = -0.5f*grid_width + k*cell_size;
                    glm::vec3 p(x, y, z);
#ifdef PERTURB
                    p += cell_size*glm::vec3(rand_float() - 0.5f, rand_float() - 0.5f, rand_float() - 0.5f);
#endif
                    grid_[i][j][k] = p;
                }
            }
        }

        // cells/tubes

        for (int i = 0; i < grid_layers; ++i)
        {
            for (int j = 0; j < grid_rows; ++j)
            {
                for (int k = 0; k < grid_cols; ++k)
                {
                    cells_.emplace_back(new sphere_geometry(grid_[i][j][k], 0.3, 3));

                    if (i < grid_layers - 1)
                        tubes_.emplace_back(new tube_geometry(grid_[i][j][k], grid_[i + 1][j][k]));

                    if (j < grid_rows - 1)
                        tubes_.emplace_back(new tube_geometry(grid_[i][j][k], grid_[i][j + 1][k]));

                    if (k < grid_cols - 1)
                        tubes_.emplace_back(new tube_geometry(grid_[i][j][k], grid_[i][j][k + 1]));
                }
            }
        }

        // lights

        for (auto &light : lights_)
        {
            light.tube = tubes_[rand_int() % tubes_.size()].get();
            light.t = rand_float();
            light.direction = rand_int() & 1;
        }
    }

    void render(const shader_program &program) const
    {
        const auto projection =
            glm::perspective(glm::radians(45.0f), static_cast<float>(window_width_) / window_height_, 0.1f, 100.f);
        const auto view_pos = glm::vec3(0, 0, -1.f - 1.f*cur_time_);
        const auto view_up = glm::vec3(cosf(0.1f*cur_time_), sinf(0.1f*cur_time_), 0);
        const auto view = glm::lookAt(view_pos, glm::vec3(0, 0, 0), view_up);

        const auto model = glm::mat4(1.0f);
        const auto mvp = projection * view * model;

        glm::mat3 model_normal(model);
        model_normal = glm::inverse(model_normal);
        model_normal = glm::transpose(model_normal);

        program.bind();
        program.set_uniform(program.uniform_location("mvp"), mvp);
        program.set_uniform(program.uniform_location("normalMatrix"), model_normal);
        program.set_uniform(program.uniform_location("modelMatrix"), model);

        std::vector<glm::vec3> light_positions;
        std::transform(lights_.begin(), lights_.end(), std::back_inserter(light_positions), [](auto &light) {
            return light.tube->pos_at(light.t);
        });
        program.set_uniform(program.uniform_location("light_positions"), light_positions);

        for (auto &cell : cells_)
            cell->render();

        for (auto &tube : tubes_)
            tube->render();
    }

    int window_width_;
    int window_height_;

    float cur_time_ = 0.f;

    shader_program program_, program_glow_;

    static constexpr auto grid_cols = 5;
    static constexpr auto grid_rows = 5;
    static constexpr auto grid_layers = 5;

    static constexpr auto cell_size = 5.f;

    static constexpr auto grid_width = cell_size * (grid_cols - 1);
    static constexpr auto grid_height = cell_size * (grid_rows - 1);
    static constexpr auto grid_depth = cell_size * (grid_layers - 1);

    std::array<std::array<std::array<glm::vec3, grid_cols>, grid_rows>, grid_layers> grid_;
    std::vector<std::unique_ptr<tube_geometry>> tubes_;
    std::vector<std::unique_ptr<sphere_geometry>> cells_;

    struct light
    {
        tube_geometry *tube;
        float t;
        bool direction;
    };
    static constexpr int num_lights = 23;
    std::array<light, num_lights> lights_;
    blur_effect blur_;
};

int main()
{
    constexpr auto window_width = 512;
    constexpr auto window_height = 512;

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

#ifdef DUMP_FRAMES
    constexpr auto total_frames = 1200;
    auto frame_num = 0;
#endif

    {
    demo d(window_width, window_height);

    while (!glfwWindowShouldClose(window))
    {
#ifdef DUMP_FRAMES
        const auto dt = 1.0f / 40;
#else
        const auto dt = 1.0f / 60;
#endif
        d.render_and_step(dt);

#ifdef DUMP_FRAMES
        char path[80];
        std::sprintf(path, "%05d.ppm", frame_num);
        dump_frame_to_file(path, window_width, window_height);

        if (++frame_num == total_frames)
            break;
#endif

        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    }

    glfwDestroyWindow(window);
    glfwTerminate();
}
