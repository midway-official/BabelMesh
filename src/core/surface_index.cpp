#include "babelmesh/cad.hpp"
#include <Eigen/Geometry>
#include <algorithm>
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace babelmesh {
void buildSurfaceIndex(SurfaceMesh &result) {
  result.bvh_nodes.clear();
  result.bvh_triangle_order.resize(result.triangles.size());
  std::iota(result.bvh_triangle_order.begin(), result.bvh_triangle_order.end(),
            0U);
  if (!result.triangles.empty()) {
    std::function<std::size_t(std::size_t, std::size_t)> buildNode;
    buildNode = [&](std::size_t begin, std::size_t end) {
      const std::size_t nodeIndex = result.bvh_nodes.size();
      result.bvh_nodes.emplace_back();
      std::array<double, 6> box{std::numeric_limits<double>::infinity(),
                                std::numeric_limits<double>::infinity(),
                                std::numeric_limits<double>::infinity(),
                                -std::numeric_limits<double>::infinity(),
                                -std::numeric_limits<double>::infinity(),
                                -std::numeric_limits<double>::infinity()};
      for (std::size_t item = begin; item < end; ++item) {
        const auto &triangle =
            result.triangles[result.bvh_triangle_order[item]];
        for (const auto vertex : triangle.vertices) {
          const auto &p = result.vertices[vertex];
          for (std::size_t axis = 0; axis < 3; ++axis) {
            box[axis] = std::min(box[axis], p[axis]);
            box[axis + 3] = std::max(box[axis + 3], p[axis]);
          }
        }
      }
      result.bvh_nodes[nodeIndex].bounds = box;
      const std::size_t count = end - begin;
      if (count <= 8) {
        result.bvh_nodes[nodeIndex].begin = begin;
        result.bvh_nodes[nodeIndex].count = count;
        return nodeIndex;
      }
      std::size_t axis = 0;
      if (box[4] - box[1] > box[3] - box[0])
        axis = 1;
      if (box[5] - box[2] > box[axis + 3] - box[axis])
        axis = 2;
      const std::size_t middle = begin + count / 2;
      auto centroid = [&](std::size_t triangleIndex) {
        const auto &triangle = result.triangles[triangleIndex];
        return (result.vertices[triangle.vertices[0]][axis] +
                result.vertices[triangle.vertices[1]][axis] +
                result.vertices[triangle.vertices[2]][axis]) /
               3.0;
      };
      std::nth_element(result.bvh_triangle_order.begin() +
                           static_cast<std::ptrdiff_t>(begin),
                       result.bvh_triangle_order.begin() +
                           static_cast<std::ptrdiff_t>(middle),
                       result.bvh_triangle_order.begin() +
                           static_cast<std::ptrdiff_t>(end),
                       [&](std::size_t lhs, std::size_t rhs) {
                         return centroid(lhs) < centroid(rhs);
                       });
      const std::size_t left = buildNode(begin, middle);
      const std::size_t right = buildNode(middle, end);
      result.bvh_nodes[nodeIndex].left = left;
      result.bvh_nodes[nodeIndex].right = right;
      return nodeIndex;
    };
    (void)buildNode(0, result.triangles.size());
  }
}

bool pointInsideSurface(const SurfaceMesh &surface,
                        const std::array<double, 3> &point) {
  if (surface.bvh_nodes.empty())
    throw std::runtime_error("surface containment requires a BVH");
  const Eigen::Vector3d origin(point[0], point[1], point[2]),
      direction =
          Eigen::Vector3d(1, .3713906763541037, .1732050807568877).normalized();
  std::vector<double> hits;
  std::vector<std::size_t> stack = {0};
  auto xyz = [&](std::size_t vertex) {
    const auto &p = surface.vertices.at(vertex);
    return Eigen::Vector3d(p[0], p[1], p[2]);
  };
  while (!stack.empty()) {
    const auto index = stack.back();
    stack.pop_back();
    const auto &node = surface.bvh_nodes.at(index);
    double near = 0, far = std::numeric_limits<double>::max();
    for (int axis = 0; axis < 3; ++axis) {
      double a = (node.bounds[static_cast<std::size_t>(axis)] - origin[axis]) /
                 direction[axis];
      double b =
          (node.bounds[static_cast<std::size_t>(axis + 3)] - origin[axis]) /
          direction[axis];
      if (a > b)
        std::swap(a, b);
      near = std::max(near, a);
      far = std::min(far, b);
    }
    if (near > far)
      continue;
    if (!node.count) {
      stack.push_back(node.left);
      stack.push_back(node.right);
      continue;
    }
    for (std::size_t i = node.begin; i < node.begin + node.count; ++i) {
      const auto &triangle =
          surface.triangles.at(surface.bvh_triangle_order.at(i));
      const Eigen::Vector3d a = xyz(triangle.vertices[0]),
                            ab = xyz(triangle.vertices[1]) - a,
                            ac = xyz(triangle.vertices[2]) - a;
      const auto cross = direction.cross(ac);
      const double determinant = ab.dot(cross);
      if (std::abs(determinant) < 1e-15)
        continue;
      const auto relative = origin - a;
      const double u = relative.dot(cross) / determinant;
      if (u < -1e-12 || u > 1 + 1e-12)
        continue;
      const auto q = relative.cross(ab);
      const double v = direction.dot(q) / determinant;
      if (v < -1e-12 || u + v > 1 + 1e-12)
        continue;
      const double distance = ac.dot(q) / determinant;
      if (distance > 1e-10)
        hits.push_back(distance);
    }
  }
  std::sort(hits.begin(), hits.end());
  std::size_t crossings = 0;
  double previous = -1;
  for (double hit : hits)
    if (previous < 0 || hit - previous > 1e-9 * std::max(1.0, std::abs(hit))) {
      ++crossings;
      previous = hit;
    }
  return crossings % 2 == 1;
}
} // namespace babelmesh
