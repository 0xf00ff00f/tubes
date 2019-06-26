#version 450 core

uniform vec3 global_light;
uniform vec3 light_positions[17];

in vec3 world_position;
in vec3 world_normal;

out vec4 frag_color;

const vec3 ambient = vec3(0.0125);
const vec3 global_light_color = vec3(0.125);
const vec3 light_color = vec3(1.0, 0.15, 0.15);

void main(void)
{
    float l = 20.0;
    for (int i = 0; i < light_positions.length(); ++i)
        l = min(l, length(light_positions[i] - world_position));
    float x = 1.0 - smoothstep(1.5, 0.25, l);

    vec3 c = max(dot(world_normal, normalize(global_light - world_position)), 0.0)*global_light_color;
    for (int i = 0; i < light_positions.length(); ++i)
    {
        vec3 d = light_positions[i] - world_position;
        float l = length(d);
        c += max(dot(world_normal, normalize(d)), 0.0)*0.05*light_color;
    }
    c += ambient;

    frag_color = vec4(mix(light_color, c, x), 1.0);
}
