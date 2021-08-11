//
// Author     : Whit Armstrong (warmstrong@anl.gov)
//
#include <XML/Helper.h>
#include "TMath.h"
#include "TString.h"
#include "DDRec/Surface.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "GeometryHelpers.h"
#include "Math/Vector3D.h"
#include "Math/AxisAngle.h"
#include "Math/VectorUtil.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::rec;

using Placements = vector<PlacedVolume>;

static Ref_t createDetector(Detector& description, xml::Handle_t e, SensitiveDetector sens){

  xml_det_t      x_det    = e;
  Material       air      = description.material("AirOptical");
  Material       vacuum   = description.vacuum();
  string         det_name = x_det.nameStr();
  string  mrichInfo = x_det.nameStr();
  //xml::Component pos      = x_det.position();
  DetElement     sdet(det_name, x_det.id());
  Assembly       assembly(det_name);
  sens.setType("photoncounter");
  OpticalSurfaceManager surfMgr = description.surfaceManager();

  // read module positions
  std::string line;
  std::ifstream input("/eic/u/sar15/eic/athena/params/projMRICHpos2.txt");
  if (!input) {
    std::cout << "MRICH module Error: file \"" << mrichInfo << "\" cannot be read." << std::endl;
  }
  int itr=0, tow = -1;
  const int nmod=124;
  double x0 = 0., y0 = 0.;
  double xx[nmod], yy[nmod];
  while (input >> tow >> x0 >> y0 )
    {
      if (!input.good()) break;
      if(itr>=nmod) break;
      xx[itr] = 0.1*x0;
      yy[itr] = 0.1*y0;
      printf("**1** %d %8.2f %8.2f\n",tow,xx[itr],yy[itr]);
      itr++;
    }

  /*
  double xx=0, yy=0, zz=0, thx=, thy=0;
  while (std::getline(input, line).good()) {
    std::istringstream iss(line);
    iss >> r >> z >> br >> bz;
    GetIndices(r, z, ir, iz, dr, dz);
    if (ir < 0 || iz < 0) {
      std::cout << "FieldMapBrBz Warning: coordinates out of range ("
		<< r << ", " << z << "), skipped it." << std::endl;
    } else {
      Bvals[ir][iz] = {br*scale, bz*scale};
      // ROOT::Math::XYZPoint p(r, 0, z);
      // std::cout << p << " -> " << trans*p << std::endl;
      // std::cout << ir << ", " << iz << ", " << br << ", " << bz << std::endl;
    }
  }
  */
  //====================================



  bool projective = getAttrOrDefault(x_det, _Unicode(projective), false);

  PlacedVolume pv;

  map<string, Volume>     modules;
  map<string, Placements> sensitives;
  map<string, Volume>     module_assemblies;
  std::map<std::string,DetElement> module_assembly_delements;

  int                     n_sensor = 1;

  xml::Component dims   = x_det.dimensions();
  auto           rmin   = dims.rmin();
  auto           rmax   = dims.rmax();
  auto           length = dims.length();
  auto           zmin   = dims.zmin();

  cout<<"******************\t"<<rmin<<"\t"<<rmax<<"\t"<<length<<"\t"<<zmin<<endl;

  // expect only one module (for now)
  xml_comp_t x_mod = x_det.child(_U(module));
  string     mod_name            = x_mod.nameStr();
  double     mod_width           = getAttrOrDefault(x_mod, _U(width), 130.0 * mm);
  double     mod_height          = getAttrOrDefault(x_mod, _U(height), 130.0 * mm);
  double     mod_length          = getAttrOrDefault(x_mod, _U(length), 130.0 * mm);

  // various components
  xml_comp_t x_frame    = x_mod.child(_Unicode(frame));
  xml_comp_t x_aerogel  = x_mod.child(_Unicode(aerogel));
  xml_comp_t x_lens     = x_mod.child(_Unicode(lens));
  xml_comp_t x_mirror   = x_mod.child(_Unicode(mirror));
  xml_comp_t x_photodet = x_mod.child(_Unicode(photodet));

  // module
  Box    m_solid(mod_width / 2.0, mod_height / 2.0, mod_length / 2.0);
  Volume m_volume(mod_name, m_solid, air);
  m_volume.setVisAttributes(description.visAttributes(x_mod.visStr()));
  DetElement mod_de( mod_name + std::string("_mod_") + std::to_string(1), 1);

  // todo module frame
  double     frame_thickness  = getAttrOrDefault(x_frame, _U(thickness), 2.0 * mm);

  // aerogel box
  xml_comp_t x_aerogel_frame = x_aerogel.child(_Unicode(frame));
  double     aerogel_width   = getAttrOrDefault(x_aerogel, _U(width), 130.0 * mm);
  double     aerogel_length  = getAttrOrDefault(x_aerogel, _U(length), 130.0 * mm);
  double     foam_thickness  = getAttrOrDefault(x_aerogel_frame, _U(thickness), 2.0 * mm);
  Material   foam_mat        = description.material(x_aerogel_frame.materialStr());
  Material   aerogel_mat     = description.material(x_aerogel.materialStr());
  auto       aerogel_vis     = getAttrOrDefault<std::string>(x_aerogel, _U(vis), std::string("InvisibleWithDaughters"));
  auto       foam_vis        = getAttrOrDefault<std::string>(x_aerogel_frame, _U(vis), std::string("RedVis"));

  // aerogel foam frame
  Box foam_box(aerogel_width / 2.0 + foam_thickness, aerogel_width / 2.0 + foam_thickness, (aerogel_length + foam_thickness) / 2.0);
  Box aerogel_sub_box(aerogel_width / 2.0, aerogel_width / 2.0, (aerogel_length + foam_thickness) / 2.0);
  SubtractionSolid foam_frame_solid(foam_box, aerogel_sub_box, Position(0, 0, foam_thickness));
  Volume           foam_vol(mod_name+"_aerogel_frame", foam_frame_solid, foam_mat);
  foam_vol.setVisAttributes(description.visAttributes(foam_vis));
  double foam_frame_zpos = -mod_length / 2.0 + frame_thickness + (aerogel_length + foam_thickness) / 2.0;
  m_volume.placeVolume(foam_vol,Position(0,0,foam_frame_zpos));

  // aerogel
  Box              aerogel_box(aerogel_width / 2.0, aerogel_width / 2.0, (aerogel_length) / 2.0);
  Volume           aerogel_vol(mod_name+"_aerogel", aerogel_box, aerogel_mat);
  aerogel_vol.setVisAttributes(description.visAttributes(aerogel_vis));
  double aerogel_zpos = foam_frame_zpos + foam_thickness / 2.0;
  pv = m_volume.placeVolume(aerogel_vol,Position(0,0,aerogel_zpos));
  DetElement aerogel_de(mod_de, mod_name + std::string("_aerogel_de") + std::to_string(1), 1);
  aerogel_de.setPlacement(pv);

  auto aerogel_surf = surfMgr.opticalSurface(dd4hep::getAttrOrDefault(x_aerogel, _Unicode(surface), "MRICH_AerogelOpticalSurface"));
  SkinSurface skin0(description, aerogel_de, Form("MRICH_aerogel_skin_surface_%d", 1), aerogel_surf, aerogel_vol);
  skin0.isValid();

  // Fresnel Lens
  //  - The lens has a constant groove pitch (delta r) as opposed to fixing the groove height.
  //  - The lens area outside of the effective diamtere is flat.
  //  - The grooves are not curved, rather they are polycone shaped, ie a flat approximating the curvature.
  auto   lens_vis       = getAttrOrDefault<std::string>(x_lens, _U(vis), std::string("AnlBlue"));
  double groove_pitch   = getAttrOrDefault(x_lens, _Unicode(pitch), 0.2 * mm);// 0.5 * mm);
  double lens_f         = getAttrOrDefault(x_lens, _Unicode(focal_length), 6.0*2.54*cm);
  double eff_diameter   = getAttrOrDefault(x_lens, _Unicode(effective_diameter), 152.4 * mm);
  //double eff_diameter   = getAttrOrDefault(x_lens, _Unicode(effective_diameter), 130.0 * mm);
  double lens_width     = getAttrOrDefault(x_lens, _Unicode(width), 6.7*2.54*cm);
  double center_thickness = getAttrOrDefault(x_lens, _U(thickness), 0.068 * 2.54 * cm);//2.0 * mm);

  double n_acrylic        = 1.49;
  double lens_curvature   = 1.0 / (lens_f*(n_acrylic - 1.0)); //confirmed
  double full_ring_rmax   = std::min(eff_diameter / 2.0, lens_width/2.0);

  double N_grooves        = std::ceil((full_ring_rmax) / groove_pitch);
  double groove_last_rmin = (N_grooves - 1) * groove_pitch;
  double groove_last_rmax = N_grooves * groove_pitch;

  auto   groove_sagitta = [&](double r) { return lens_curvature * std::pow(r, 2) / (1.0 + 1.0); };
  double lens_thickness = groove_sagitta(groove_last_rmax) - groove_sagitta(groove_last_rmin) + center_thickness;

  Material         lens_mat = description.material(x_lens.materialStr());
  Box              lens_box(lens_width / 2.0, lens_width / 2.0, (center_thickness) / 2.0);
  SubtractionSolid flat_lens(lens_box, Tube(0.0, full_ring_rmax, 2 * center_thickness));

  Assembly lens_vol(mod_name + "_lens");
  Volume   flatpart_lens_vol( "flatpart_lens", flat_lens, lens_mat);
  lens_vol.placeVolume(flatpart_lens_vol);//,Position(0,0,lens_zpos));

  Solid  fresnel_lens_solid;

  int    i_groove           = 0;
  double groove_rmax        = groove_pitch;
  double groove_rmin        = 0;

  cout<<"full_ring_rmax: \t"<<full_ring_rmax<<endl;

  while ( groove_rmax <= full_ring_rmax ) {
    double   dZ = groove_sagitta(groove_rmax) - groove_sagitta(groove_rmin);
    //std::cout << " dZ = " << dZ << ",  lens_thickness = " << lens_thickness << "\n";
    //std::cout  << "groove_rmin = " << groove_rmin << "\n";
    //std::cout  << "groove_rmax = " << groove_rmax << "\n";
    Polycone groove_solid(0, 2.0 * M_PI,
                          {groove_rmin, groove_rmin, groove_rmin},
                          {groove_rmax, groove_rmax, groove_rmin},
                          {-lens_thickness/2.0, lens_thickness/2.0-dZ, lens_thickness/2.0});
    Volume   lens_groove_vol("lens_groove_" + std::to_string(i_groove), groove_solid, lens_mat);
    lens_vol.placeVolume(lens_groove_vol); //,Position(0,0,lens_zpos));
    //Volume   groove_vol(groove_solid, lens_mat, par->name.c_str(), 0, 0, 0);
    //new G4PVPlacement(0, par->pos, Groove_log[i], par->name.c_str(), motherLV, false, 0, OverlapCheck());
    //phi1 = phi1 + halfpi;  //g4 pre-defined: halfpi=pi/2
    //Tube       sub_cylinder(r0, r1, 3*eff_diameter);
    //IntersectionSolid groove_solid(lens_box,lens_sphere, Position(0,0,-eff_diameter/2.0 + lens_thickness/2.0+(t-lens_t)/2.0 ));
    //IntersectionSolid lens_ring(groove_solid, sub_cylinder);
    //if (i_groove == 0) {
    //  fresnel_lens_solid = groove_solid;
    //} else {
    //  fresnel_lens_solid = UnionSolid(fresnel_lens_solid, groove_solid);
    //}
    //r0 = r1;
    //if(i_groove > 3) { 
    //  SubtractionSolid flat_lens(lens_box,Tube(0.0, r0, 3*eff_diameter));
    //  fresnel_lens_solid = UnionSolid(fresnel_lens_solid, flat_lens);
    //  break; // temporary
    //}
    i_groove++;
    groove_rmin = (i_groove  )*groove_pitch;
    groove_rmax = (i_groove+1)*groove_pitch;
  }

  cout<<"N_grooves :"<<i_groove<<endl;
  //fresnel_lens_solid = UnionSolid(fresnel_lens_solid, flat_lens);
  //Volume   lens_vol(mod_name + "_lens", fresnel_lens_solid, lens_mat);

  lens_vol.setVisAttributes(description.visAttributes(lens_vis));
  double lens_zpos = aerogel_zpos +aerogel_length/ 2.0 + foam_thickness + lens_thickness/2.0;
  pv = m_volume.placeVolume(lens_vol,Position(0,0,lens_zpos));
  DetElement lens_de(mod_de, mod_name + std::string("_lens_de") + std::to_string(1), 1);
  lens_de.setPlacement(pv);

  auto surf = surfMgr.opticalSurface(dd4hep::getAttrOrDefault(x_lens, _Unicode(surface), "MRICH_LensOpticalSurface"));
  SkinSurface skin(description, lens_de, Form("MRichFresnelLens_skin_surface_%d", 1), surf, lens_vol);
  skin.isValid();

  // mirror
  auto   mirror_vis      = getAttrOrDefault<std::string>(x_mirror, _U(vis), std::string("AnlGray"));
  double mirror_x1      = getAttrOrDefault(x_mirror, _U(x1), 100.0 * mm);
  double mirror_x2      = getAttrOrDefault(x_mirror, _U(x2), 80.0 * mm);
  double mirror_length  = getAttrOrDefault(x_mirror, _U(length), 130.0 * mm);
  double mirror_thickness  = getAttrOrDefault(x_mirror, _U(thickness), 2.0 * mm);
  double outer_x1 = (mirror_x1+mirror_thickness)/2.0;
  double outer_x2 = (mirror_x2+mirror_thickness)/2.0;
  Trd2   outer_mirror_trd(outer_x1, outer_x2, outer_x1,  outer_x2, mirror_length/2.0);
  Trd2   inner_mirror_trd(mirror_x1 / 2.0,  mirror_x2 / 2.0, mirror_x1 / 2.0,mirror_x2 / 2.0, mirror_length/2.0+0.1*mm);
  SubtractionSolid mirror_solid(outer_mirror_trd, inner_mirror_trd);
  Material   mirror_mat        = description.material(x_mirror.materialStr());
  Volume           mirror_vol(mod_name+"_mirror", mirror_solid, mirror_mat);
  double mirror_zpos = lens_zpos + lens_thickness/2.0 + foam_thickness + mirror_length/2.0;
  pv = m_volume.placeVolume(mirror_vol,Position(0,0,mirror_zpos));
  DetElement mirror_de(mod_de, mod_name + std::string("_mirror_de") + std::to_string(1), 1);
  mirror_de.setPlacement(pv);

  auto mirror_surf = surfMgr.opticalSurface(dd4hep::getAttrOrDefault(x_mirror, _Unicode(surface), "MRICH_MirrorOpticalSurface"));
  SkinSurface skin1(description, mirror_de, Form("MRICH_mirror_skin_surface_%d", 1), mirror_surf, mirror_vol);
  skin1.isValid();

  // photon detector
  xml_comp_t x_photodet_sensor  = x_photodet.child(_Unicode(sensor));
  auto   photodet_vis      = getAttrOrDefault<std::string>(x_photodet, _U(vis), std::string("AnlRed"));
  double     photodet_width     = getAttrOrDefault(x_photodet, _U(width), 130.0 * mm);
  double     photodet_thickness = getAttrOrDefault(x_photodet, _U(thickness), 2.0 * mm);
  double     sensor_thickness   = getAttrOrDefault(x_photodet_sensor, _U(thickness), 2.0 * mm);
  Material   photodet_mat       = description.material(x_photodet.materialStr());
  Material   sensor_mat         = description.material(x_photodet_sensor.materialStr());
  int     sensor_nx   = getAttrOrDefault(x_photodet_sensor, _Unicode(nx), 2);
  int     sensor_ny   = getAttrOrDefault(x_photodet_sensor, _Unicode(ny), 2);

  Box    window_box(photodet_width/2.0,photodet_width/2.0,photodet_thickness/2.0);
  Volume           window_vol(mod_name+"_window", window_box, photodet_mat);
  window_vol.setSensitiveDetector(sens);
  double window_zpos = mirror_zpos + mirror_length/2.0+photodet_thickness/2.0;
  pv = m_volume.placeVolume(window_vol,Position(0,0,window_zpos));
  DetElement   comp_de(mod_de, std::string("mod_sensor_de_") + std::to_string(1) ,  1);
  comp_de.setPlacement(pv);
  pv.addPhysVolID("sensor", n_sensor);

  //for (size_t ic = 0; ic < sensVols.size(); ++ic) {
  //  PlacedVolume sens_pv = sensVols[ic];
  //  DetElement   comp_de(mod_de, std::string("de_") + sens_pv.volume().name(), ic + 1);
  //  comp_de.setPlacement(sens_pv);
  //  // Acts::ActsExtension* sensorExtension = new Acts::ActsExtension();
  //  //// sensorExtension->addType("sensor", "detector");
  //  // comp_de.addExtension<Acts::ActsExtension>(sensorExtension);
  //  //// comp_de.setAttributes(description, sens_pv.volume(), x_layer.regionStr(),
  //  //// x_layer.limitsStr(),
  //  ////                      xml_det_t(xmleles[m_nam]).visStr());
  //}
  //DetElement window_de(sdet, mod_name + std::string("_window_de") + std::to_string(1), 1);
  //window_de.setPlacement(pv);

  window_vol.setSensitiveDetector(sens);
  sensitives[mod_name].push_back(pv);
  ++n_sensor;
  modules[mod_name]                   = m_volume;
  module_assembly_delements[mod_name] = mod_de;
  // end module

  int i_mod = 1;
  // detector envelope
  Tube   envShape(rmin, rmax, length / 2., 0., 2 * M_PI);
  Volume envVol("MRICH_Envelope", envShape, air);
  envVol.setVisAttributes(description.visAttributes(x_det.visStr()));

  // place modules in the sectors (disk)
  auto points = athena::geo::fillSquares({0., 0.}, mod_width, rmin, rmax);

  // mod_name = ...
  Placements& sensVols = sensitives[mod_name];
  auto        mod_v    = modules[mod_name];
  // determine module direction, always facing z = 0
  double roty = dims.zmin() < 0. ? -M_PI : 0 ;

  /*
  int    imod = 1;
  for (auto& p : points) {


  }
  */
  //*
  //cout<<"*********\t"<<points<<"\t*********"<endl;
  double xl[]={0.1,0.3,0.8};
  int  ic = 0,  imod = 1;
  double scale = 1.0;
  for(int imod=0; imod<nmod; imod++){
  //for (auto& p : points) {
    if(ic<8||ic>107) continue;
    //if(imod<44){
    //ROOT::Math::XYZVector x_location(xx[ic], yy[ic], zz[ic]);
     //} else 

     /*
    ROOT::Math::XYZVector x_location(p.x(), p.y(), zmin+std::signbit(zmin)*mod_length/2.0);
    ROOT::Math::XYZVector z_dir(0, 0, 1);
    ROOT::Math::XYZVector x_dir(1, 0, 0);
    ROOT::Math::XYZVector rot_axis = x_location.Cross(z_dir);
    double rot_angle = ROOT::Math::VectorUtil::Angle(z_dir,x_location);
    ROOT::Math::AxisAngle proj_rot(rot_axis,-1.0*rot_angle);
    ROOT::Math::AxisAngle grid_fix_rot(x_location,0.0*rot_angle);
    auto new_x_dir = grid_fix_rot*x_dir;
    // operations are inversely ordered
    Transform3D tr = Translation3D(p.x(), p.y(), 0.) // move to position
                     * RotationY(roty);              // facing z = 0.

    */
    if(ic<8) scale = 1.062;
    else if(ic>=8 && ic<20) scale = 1.065;
    else if(ic>=20 && ic<24) scale = 1.067;
    else if(ic>=24 && ic<44) scale = 1.071;
    else if(ic>=44 && ic<48) scale = 1.075;
    else if(ic>=48 && ic<76) scale = 1.084;
    else if(ic>=76 && ic<80) scale = 1.088;
    else if(ic>=80 && ic<108) scale = 1.092;
    else scale = 1.11;

    xx[ic]*=scale;
    yy[ic]*=scale;

    double zz = -175.0;
    double rotAngX = atan(yy[ic]/zz);//*180./3.14159;
    double rotAngY = -1.*atan(xx[ic]/zz);//*180./3.14159;

    Transform3D tr = Translation3D(xx[ic], yy[ic], zz)*RotationX(rotAngX)*RotationY(rotAngY);
    /*
    if(projective) {
      tr = Translation3D(p.x(), p.y(), 0.)              // move to position
           * grid_fix_rot  // keep the modules oriented vertially
           * proj_rot      // projective rotation
           * RotationY(roty); // facing z = 0.
    }
    */
    // mod placement
    pv = envVol.placeVolume(mod_v, tr);
    pv.addPhysVolID("module", i_mod);

    auto mod_det_element  = module_assembly_delements[mod_name].clone(mod_name + "__" + std::to_string(i_mod));
    mod_det_element.setPlacement(pv);
    sdet.add(mod_det_element);
    //cout<<"........MOD.........\t"<<i_mod<<"\t"<<xx[ic]<<"\t"<<yy[ic]<<"\t"<<zz[ic]<<endl;//"\t"<<p.x()<<"\t"<<p.y() <<endl;
    i_mod++;
    ic++;
  }
  cout<<i_mod<<endl;
  //*/

  // place envelope
  Volume       motherVol = description.pickMotherVolume(sdet);
  PlacedVolume envPV     = motherVol.placeVolume(envVol, Position(0, 0, zmin));
  envPV.addPhysVolID("system", x_det.id());
  sdet.setPlacement(envPV);
  //cout<<"*********************"<<zmin <<"*********************"<<endl;
  return sdet;
}

//void addModules(Volume &mother, xml::DetElement &detElem, Detector &description, SensitiveDetector &sens)
//{
//    xml::Component dims = detElem.dimensions();
//    xml::Component mods = detElem.child(_Unicode(modules));
//
//    auto rmin = dims.rmin();
//    auto rmax = dims.rmax();
//
//    auto mThick = mods.attr<double>(_Unicode(thickness));
//    auto mWidth = mods.attr<double>(_Unicode(width));
//    auto mGap = mods.attr<double>(_Unicode(gap));
//
//    auto modMat = description.material(mods.materialStr());
//    auto gasMat = description.material("AirOptical");
//
//    // single module
//    Box mShape(mWidth/2., mWidth/2., mThick/2. - 0.1*mm);
//    Volume mVol("ce_MRICH_mod_Solid", mShape, modMat);
//
//    // a thin gas layer to detect optical photons
//    Box modShape(mWidth/2., mWidth/2., mThick/2.);
//    Volume modVol("ce_MRICH_mod_Solid_v", modShape, gasMat);
//    // thin gas layer is on top (+z) of the material
//    modVol.placeVolume(mVol, Position(0., 0., -0.1*mm));
//
//    modVol.setVisAttributes(description.visAttributes(mods.visStr()));
//    sens.setType("photoncounter");
//    modVol.setSensitiveDetector(sens);
//
//    // place modules in the sectors (disk)
//    auto points = ref::utils::fillSquares({0., 0.}, mWidth + mGap, rmin - mGap, rmax + mGap);
//
//    // determine module direction, always facing z = 0
//    double roty = dims.z() > 0. ? M_PI/2. : -M_PI/2.;
//    int imod = 1;
//    for (auto &p : points) {
//        // operations are inversely ordered
//        Transform3D tr = Translation3D(p.x(), p.y(), 0.)        // move to position
//                       * RotationY(roty);                       // facing z = 0.
//        auto modPV = mother.placeVolume(modVol, tr);
//        modPV.addPhysVolID("sector", 1).addPhysVolID("module", imod ++);
//    }
//}

// clang-format off
DECLARE_DETELEMENT(athena_MRICH, createDetector)

