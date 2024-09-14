#include <iostream>
#include <vector>
#include "Common.h"
#include "Scene.h"
#include "Camera.h"
#include "Material.h"
#include "Hittable.h"
#include "Utils/lodepng.h"

const int kMaxTraceDepth = 5;

Color TraceRay(const Ray& ray,
	const std::vector<LightSource>& light_sources,
	const Hittable& scene,
	int trace_depth);

Color Shade(const std::vector<LightSource>& light_sources,
	const Hittable& hittable_collection,
	const HitRecord& hit_record,
	int trace_depth) {

	Color color(0.f, 0.f, 0.f);
	Material m = hit_record.material;

	//Ambient term
	color = m.k_a * m.ambient;

	//For each light source, add Diffuse and specular terms due to each light source to color
	for (int i = 0;i < light_sources.size();i++) {
		LightSource light = light_sources[i];
		Ray shadow_ray(hit_record.position, glm::normalize(light.position - hit_record.position));
		HitRecord hr;

		float l_dot_n = glm::dot(hit_record.normal, shadow_ray.d);
		if (l_dot_n >= 0) {
			if (!hittable_collection.Hit(shadow_ray, &hr)) { // if shadow_ray is NOT blocked
				Vec r = glm::normalize(2 * l_dot_n * hit_record.normal - shadow_ray.d); //direction of the reflection of shadow ray (normalized)
				float r_dot_v = glm::dot(r, -hit_record.in_direction);
				if (r_dot_v < 1e-6) r_dot_v = 0;
				color.r += light.intensity * (m.k_d * m.diffuse.r * l_dot_n + m.k_s * m.specular.r * pow(r_dot_v, m.sh));
				color.g += light.intensity * (m.k_d * m.diffuse.g * l_dot_n + m.k_s * m.specular.g * pow(r_dot_v, m.sh));
				color.b += light.intensity * (m.k_d * m.diffuse.b * l_dot_n + m.k_s * m.specular.b * pow(r_dot_v, m.sh));
			}
		}

	}

	//Reflected ray contribution
	if (trace_depth < kMaxTraceDepth) {
		if (m.k_s > 0) {
			Ray reflected_ray(hit_record.position, hit_record.reflection);
			Color r_color = TraceRay(reflected_ray, light_sources, hittable_collection, trace_depth + 1);
			color += r_color * m.k_s; //scale r_color by reflectance and add to color;
		}
	}

	//clamp color to 1.0
	color.r = color.r > 1.0 ? 1. : color.r;
	color.g = color.g > 1.0 ? 1. : color.g;
	color.b = color.b > 1.0 ? 1. : color.b;

	return color;
}

Color TraceRay(const Ray& ray,
	const std::vector<LightSource>& light_sources,
	const Hittable& hittable_collection,
	int trace_depth) {
	HitRecord record;
	Color color(0.0f, 0.0f, 0.0f);
	if (hittable_collection.Hit(ray, &record)) {
		return Shade(light_sources, hittable_collection, record, trace_depth);
	}
	return color; //if no hit, return background color (0,0,0)
}

int main() {
	// TODO: Set your workdir (absolute path) here.
    const std::string work_dir("D:/COMP3271_PA3/");

	// Construct scene
    Scene scene(work_dir, "scene/sphere.toml");
	const Camera& camera = scene.camera_;
	int width = camera.width_;
	int height = camera.height_;

	std::vector<unsigned char> image(width * height * 4, 0);

	float progress = 0.f;

	// Traverse all pixels
	for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; y++) {
			Color color(0.f, 0.f, 0.f);
			int count = 0;
			for (float bias_x = 0.25f; bias_x < 1.f; bias_x += .5f) {
				for (float bias_y = 0.25f; bias_y < 1.f; bias_y += .5f) {
					Ray ray = camera.RayAt(float(x) + bias_x, float(y) + bias_y);
					color += TraceRay(ray, scene.light_sources_, scene.hittable_collection_, 1);
					count++;
				}
			}
			color /= float(count);
			int idx = 4 * ((height - y - 1) * width + x);
			for (int i = 0; i < 3; i++) {
				image[idx + i] = (uint8_t)(glm::min(color[i], 1.f - 1e-5f) * 256.f);
			}
			image[idx + 3] = 255;

			float curr_progress = float(x * height + y) / float(height * width);
			if (curr_progress > progress + 0.05f) {
				progress += 0.05f;
				std::cout << "Progress: " << progress << std::endl;
			}
		}
	}

	// Save result as png file
	std::vector<unsigned char> png;
	unsigned error = lodepng::encode(png, image, width, height);
	lodepng::save_file(png, work_dir + "output.png");
}
