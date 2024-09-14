#include "Hittable.h"



// Sphere
bool Sphere::Hit(const Ray& ray, HitRecord* hit_record) const {

	if (glm::distance(ray.o, this->o_) < this->r_ + 1e-5f) { //if ray.o is inside the sphere, return false
		return false;
	}

	//Consider a sphere centered as the origin for calculation
	Point fake_o = ray.o - this->o_;

	//Quadratic equation
	float a = 1; //A = d.d = 1
	float b = 2 * glm::dot(fake_o, ray.d); // B = 2o.d
	float r = this->r_; 
	float c = glm::dot(fake_o, fake_o) - r * r; // C = o.o -r^2
	float d = b * b - 4 * a * c; //Discriminant = b^2-4ac

	if (d >= 0.f) {
		float t1 = (-b + sqrt(d)) / (2 * a); //solve for quadratic equation  - root 1
		float t2 = (-b - sqrt(d)) / (2 * a); //solve for quadratic equation  - root 2
		float t = t1 > t2 ? t2 : t1; //get the smaller root
		if (t < 0)return false; // return false if the "intersection" is in the opposite direction of the ray
		hit_record->position = ray.At(t); //intersected position
		hit_record->normal = glm::normalize(hit_record->position - this->o_); //intersected surface normal vector = normalize[vector from sphere center to intersection position]
		hit_record->distance = glm::distance(hit_record->position, ray.o); //distance from ray’s origin to intersection;
		hit_record->in_direction = ray.d;
		hit_record->reflection = glm::normalize(ray.d - 2 * glm::dot(ray.d, hit_record->normal) * hit_record->normal);//reflected direction(normalized)
		hit_record->material = material_;

		return true;
	}
	return false;
}

// Quadric
bool Quadric::Hit(const Ray& ray, HitRecord* hit_record) const {
	glm::vec4 rayO(ray.o[0], ray.o[1], ray.o[2], 1);
	glm::vec4 rayD(ray.d[0], ray.d[1], ray.d[2], 0);

	// We can calculate (x^T A x) using glm::dot(x,A*x)

	//Quadratic equation
	float a = glm::dot(rayD, this->A_ * rayD); //A = (D^T A D)
	float b = 2 * glm::dot(rayO, this->A_ * rayD); // B = 2(O^T A D)
	float c = glm::dot(rayO, this->A_ * rayO); // C = (O^T A O)
	float d = b * b - 4 * a * c; //Discriminant = b^2-4ac

	if (d >= 0) {
		float t1 = (-b + sqrt(d)) / (2 * a);  //solve for quadratic equation  - root 1
		float t2 = (-b - sqrt(d)) / (2 * a);  //solve for quadratic equation  - root 2
		float t = t1 > t2 ? t2 : t1; //get the smaller root
		if (t < 0)return false; // return false if the "intersection" is in the opposite direction of the ray
		
		hit_record->position = ray.At(t); //intersected position
		hit_record->normal = glm::normalize(Vec((this->A_ + glm::transpose(this->A_)) * glm::vec4(hit_record->position.x, hit_record->position.y, hit_record->position.z, 1)));//intersected surface normal vector = normalize[(A + A^T)x]
		hit_record->distance = glm::distance(hit_record->position, ray.o); //distance from ray’s origin to intersection;
		hit_record->in_direction = ray.d;
		hit_record->reflection = glm::normalize(ray.d - 2 * glm::dot(ray.d, hit_record->normal) * hit_record->normal);//reflected direction(normalized)
		hit_record->material = material_;
		return true;
	}
	return false;


}

// Triangle
bool Triangle::Hit(const Ray& ray, HitRecord* hit_record) const {
	
	Vec normal = glm::normalize(glm::cross(this->b_ - this->a_, this->c_ - this->a_));
	float d = glm::dot(normal, this->a_); // n.x=d
	if (glm::dot(normal, ray.d) != 0) {
		/*
		n.x=d (1)
		Ray: x = o+dt (2)
		Sub (2) into (1):
			n.(o+dt) = d
			n.o+n.dt = d
				n.dt = d-n.o
				   t = (d-n.o)/(n.d)
		*/
		float t = (d - glm::dot(normal, ray.o)) / glm::dot(normal, ray.d);
		if (t < 1e-5) return false;
		Vec p = ray.At(t);
		Vec pa = this->a_ - p; // vector from vertex a to the interestion point
		Vec pb = this->b_ - p; // vector from vertex b to the interestion point
		Vec pc = this->c_ - p; // vector from vertex c to the interestion point
		Vec x1 = glm::normalize(glm::cross(pa, pb)); 
		Vec x2 = glm::normalize(glm::cross(pb, pc)); 
		Vec x3 = glm::normalize(glm::cross(pc, pa)); 
		if (glm::dot(x1, x2) > 0 && glm::dot(x1, x3) > 0) { //if the interestion point is inside the triangle, the cross production have the same direction 
			hit_record->position = p; //intersected position

			//intersected surface normal vector
			if (phong_interpolation_) {
				float abc = glm::length(glm::cross(this->b_ - this->a_, this->c_ - this->a_));
				float a1 = glm::length(glm::cross(this->b_ - p, this->c_ - p)) / abc;
				float a2 = glm::length(glm::cross(this->c_ - p, this->a_ - p)) / abc;
				float a3 = 1.0f - a1 - a2;
				hit_record->normal = a1 * this->n_a_ + a2 * this->n_b_ + a3 * this->n_c_;
			} else {
				hit_record->normal = normal;
			}

			hit_record->distance = glm::distance(hit_record->position, ray.o); //distance from ray’s origin to intersection;
			hit_record->in_direction = ray.d;
			hit_record->reflection = glm::normalize(ray.d - 2 * glm::dot(ray.d, hit_record->normal) * hit_record->normal);//reflected direction(normalized)
			return true;
		}
	}
	return false;
}

// ---------------------------------------------------------------------------------------------
// ------------------------------ no need to change --------------------------------------------
// ---------------------------------------------------------------------------------------------

// CompleteTriangle
bool CompleteTriangle::Hit(const Ray& ray, HitRecord* hit_record) const {
	bool ret = triangle_.Hit(ray, hit_record);
	if (ret) {
		hit_record->material = material_;
	}
	return ret;
}


// Mesh
Mesh::Mesh(const std::string& file_path,
	const Material& material,
	bool phong_interpolation) :
	ply_data_(file_path), material_(material), phong_interpolation_(phong_interpolation) {
	std::vector<std::array<double, 3>> v_pos = ply_data_.getVertexPositions();
	vertices_.resize(v_pos.size());

	for (int i = 0; i < vertices_.size(); i++) {
		vertices_[i] = Point(v_pos[i][0], v_pos[i][1], v_pos[i][2]);
	}

	f_ind_ = ply_data_.getFaceIndices();

	// Calc face normals
	for (const auto& face : f_ind_) {
		Vec normal = glm::normalize(glm::cross(vertices_[face[1]] - vertices_[face[0]], vertices_[face[2]] - vertices_[face[0]]));
		face_normals_.emplace_back(normal);
	}

	// Calc vertex normals
	vertex_normals_.resize(vertices_.size(), Vec(0.f, 0.f, 0.f));
	for (int i = 0; i < f_ind_.size(); i++) {
		for (int j = 0; j < 3; j++) {
			vertex_normals_[f_ind_[i][j]] += face_normals_[i];
		}
	}
	for (auto& vertex_normal : vertex_normals_) {
		vertex_normal = glm::normalize(vertex_normal);
	}

	// Construct hittable triangles
	for (const auto& face : f_ind_) {
		triangles_.emplace_back(vertices_[face[0]], vertices_[face[1]], vertices_[face[2]],
			vertex_normals_[face[0]], vertex_normals_[face[1]], vertex_normals_[face[2]],
			phong_interpolation_);
	}

	// Calc bounding box
	Point bbox_min(1e5f, 1e5f, 1e5f);
	Point bbox_max(-1e5f, -1e5f, -1e5f);
	for (const auto& vertex : vertices_) {
		bbox_min = glm::min(bbox_min, vertex - 1e-3f);
		bbox_max = glm::max(bbox_max, vertex + 1e-3f);
	}

	// Build Octree
	tree_nodes_.emplace_back(new OctreeNode());
	tree_nodes_.front()->bbox_min = bbox_min;
	tree_nodes_.front()->bbox_max = bbox_max;

	root_ = tree_nodes_.front().get();
	for (int i = 0; i < f_ind_.size(); i++) {
		InsertFace(root_, i);
	}
}

bool Mesh::Hit(const Ray& ray, HitRecord* hit_record) const {
	const bool brute_force = false;
	if (brute_force) {
		// Naive hit algorithm
		float min_dist = 1e5f;
		for (const auto& triangle : triangles_) {
			HitRecord curr_hit_record;
			if (triangle.Hit(ray, &curr_hit_record)) {
				if (curr_hit_record.distance < min_dist) {
					*hit_record = curr_hit_record;
					min_dist = curr_hit_record.distance;
				}
			}
		}
		if (min_dist + 1.0 < 1e5f) {
			hit_record->material = material_;
			return true;
		}
		return false;
	} else {
		bool ret = OctreeHit(root_, ray, hit_record);
		if (ret) {
			hit_record->material = material_;
		}
		return ret;
	}
}

bool Mesh::IsFaceInsideBox(const std::vector<size_t>& face, const Point& bbox_min, const Point& bbox_max) const {
	for (size_t idx : face) {
		const auto& pt = vertices_[idx];
		for (int i = 0; i < 3; i++) {
			if (pt[i] < bbox_min[i] + 1e-6f) return false;
			if (pt[i] > bbox_max[i] - 1e-6f) return false;
		}
	}
	return true;
}

bool Mesh::IsRayIntersectBox(const Ray& ray, const Point& bbox_min, const Point& bbox_max) const {
	float t_min = -1e5f;
	float t_max = 1e5f;

	for (int i = 0; i < 3; i++) {
		if (glm::abs(ray.d[i]) < 1e-6f) {
			if (ray.o[i] < bbox_min[i] + 1e-6f || ray.o[i] > bbox_max[i] - 1e-6f) {
				t_min = 1e5f;
				t_max = -1e5f;
			}
        }
        else {
			if (ray.d[i] > 0.f) {
				t_min = glm::max(t_min, (bbox_min[i] - ray.o[i]) / ray.d[i]);
				t_max = glm::min(t_max, (bbox_max[i] - ray.o[i]) / ray.d[i]);
            }
            else {
				t_min = glm::max(t_min, (bbox_max[i] - ray.o[i]) / ray.d[i]);
				t_max = glm::min(t_max, (bbox_min[i] - ray.o[i]) / ray.d[i]);
			}
		}
	}

	return t_min + 1e-6f < t_max;
}

void Mesh::InsertFace(OctreeNode* u, size_t face_idx) {
	const Point& bbox_min = u->bbox_min;
	const Point& bbox_max = u->bbox_max;

	Vec bias = bbox_max - bbox_min;
	Vec half_bias = bias * 0.5f;

	bool inside_childs = false;

	for (size_t a = 0; a < 2; a++) {
		for (size_t b = 0; b < 2; b++) {
			for (size_t c = 0; c < 2; c++) {
				size_t child_idx = ((a << 2) | (b << 1) | c);
				Point curr_bbox_min = bbox_min + half_bias * Vec(float(a), float(b), float(c));
				Point curr_bbox_max = curr_bbox_min + half_bias;
				if (IsFaceInsideBox(f_ind_[face_idx], curr_bbox_min, curr_bbox_max)) {
					if (u->childs[child_idx] == nullptr) {
						tree_nodes_.emplace_back(new OctreeNode());
						OctreeNode* child = tree_nodes_.back().get();
						u->childs[child_idx] = tree_nodes_.back().get();
						child->bbox_min = curr_bbox_min;
						child->bbox_max = curr_bbox_max;
					}
					InsertFace(u->childs[child_idx], face_idx);
					inside_childs = true;
				}
			}
		}
	}

	if (!inside_childs) {
		u->face_index.push_back(face_idx);
	}
}

bool Mesh::OctreeHit(OctreeNode* u, const Ray& ray, HitRecord* hit_record) const {
	if (!IsRayIntersectBox(ray, u->bbox_min, u->bbox_max)) {
		return false;
	}
	float distance = 1e5f;
	for (const auto& face_idx : u->face_index) {
		HitRecord curr_hit_record;
		if (triangles_[face_idx].Hit(ray, &curr_hit_record)) {
			if (curr_hit_record.distance < distance) {
				distance = curr_hit_record.distance;
				*hit_record = curr_hit_record;
			}
		}
	}

	for (const auto& child : u->childs) {
		if (child == nullptr) {
			continue;
		}
		HitRecord curr_hit_record;
		if (OctreeHit(child, ray, &curr_hit_record)) {
			if (curr_hit_record.distance < distance) {
				distance = curr_hit_record.distance;
				*hit_record = curr_hit_record;
			}
		}
	}
	return distance + 1 < 1e5f;
}


// Hittable list
void HittableList::PushHittable(const Hittable& hittable) {
	hittable_list_.push_back(&hittable);
}

bool HittableList::Hit(const Ray& ray, HitRecord* hit_record) const {
	float min_dist = 1e5f;
	for (const auto& hittable : hittable_list_) {
		HitRecord curr_hit_record;
		if (hittable->Hit(ray, &curr_hit_record)) {
			if (curr_hit_record.distance < min_dist) {
				*hit_record = curr_hit_record;
				min_dist = curr_hit_record.distance;
			}
		}
	}
	return min_dist + 1.0 < 1e4f;
}