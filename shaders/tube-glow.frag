#version 450 core

uniform vec3 light_positions[23];

in vec3 world_position;
in vec3 world_normal;

out vec4 frag_color;

const vec3 light_color = vec3(1.0, 0.5, 0.5);

void main(void)
{
    float l = 20.0;
    for (int i = 0; i < light_positions.length(); ++i)
        l = min(l, length(light_positions[i] - world_position));
    float x = smoothstep(1.5, 0.25, l);
    frag_color = vec4(light_color*x, 1.0);
}
