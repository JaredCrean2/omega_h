#include "Omega_h_mesh.hpp"
#include "Omega_h_metric.hpp"
#include "Omega_h_adapt.hpp"
#include "Omega_h_file.hpp"
#include "Omega_h_array_ops.hpp"
#include "Omega_h_quality.hpp"

using Real = Omega_h::Real;

void refine_metric(Omega_h::Write<Real>& metric, Real factor)
{
  assert(metric.size() % 6 == 0);
  int nverts = metric.size()/6;
  for (int i=0; i < nverts; ++i)
    for (int j=0; j < 6; ++j)
      metric[6*i + j] *= factor*factor;
}

void run_refinement(Omega_h::Mesh* mesh, const std::string& outname_prefix, double refine_factor)
{
  std::cout << "number of verts in initial mesh = " << mesh->nverts() << std::endl;
  std::cout << "number of elements in initial mesh = " << mesh->nelems() << std::endl;

  std::string metric_field_name = "new_metric";
  Omega_h::Reals orig_metric = Omega_h::get_implied_metrics(mesh);
  Omega_h::Write<Real> target_metric = Omega_h::deep_copy(orig_metric, "target_metric");
  refine_metric(target_metric, refine_factor);
  mesh->add_tag(Omega_h::VERT, metric_field_name, 6, Omega_h::Read<Real>(target_metric));

  Omega_h::vtk::write_vtu(outname_prefix + "_initial_mesh.vtu", mesh);


  int num_initial_verts = mesh->nverts();
  Omega_h::AdaptOpts opts(mesh);
  opts.min_quality_allowed = 0.15;  // this is needed for the 3X refine case
  opts.xfer_opts.type_map[metric_field_name] = Omega_h_Transfer::OMEGA_H_METRIC;

  Omega_h::grade_fix_adapt(mesh, opts, target_metric, true);

  std::cout << "number of verts in new mesh = " << mesh->nverts() << std::endl;
  std::cout << "number of elements in new mesh = " << mesh->nelems() << std::endl;

  Omega_h::Read<Real> target_metric2 = mesh->get_array<Real>(0, metric_field_name);
  std::cout << "min quality = " << Omega_h::get_min(mesh->comm(), Omega_h::measure_qualities(mesh, target_metric2)) << std::endl;

  Omega_h::vtk::write_vtu(outname_prefix + "_adapted_mesh.vtu", mesh);
}

int main(int argc, char** argv)
{
  std::string fname = "stv.osh";
  
  auto lib = Omega_h::Library(&argc, &argv);
  Omega_h::Mesh mesh(&lib);
  Omega_h::binary::read(fname, lib.world(), &mesh);  
  
  run_refinement(&mesh, "stv_coarsen2", 0.5);
  
  return 0;
}