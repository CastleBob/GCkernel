//! Independent construction/transport/display acceptance for full-turn polygonal revolve.
//! No shared worker is rebuilt by these tests. Dimensions and expectations are analytic.
use std::collections::HashMap;
use std::f64::consts::PI;
use gc_kernel::geom::{Surface, Vec3};
use gc_kernel::neutral::{Body, Reader, Profile, write_body};
use gc_kernel::revolve::{Section, revolve};
use gc_kernel::{mass, protocol, tessellate, topology};

const TOL: f64 = 1e-7;
fn section(points: Vec<(f64,f64)>) -> Section {
    Section { p: Vec3::new(0.,0.,0.), a: Vec3::new(0.,1.,0.), x: Vec3::new(1.,0.,0.), points }
}
fn close(a: f64, b: f64) {
    assert!(a.is_finite() && (a-b).abs() <= 1e-8*b.abs().max(1.), "{a} vs analytic {b}");
}
fn frustum(a: f64, b: f64, h: f64) -> f64 { PI*h*(a*a+a*b+b*b)/3. }

fn checked_mesh(body: &Body, chord: f64) {
    topology::check(body,TOL).unwrap();
    let wire = write_body(body);
    let mut reader = Reader::with_tolerance(&wire,TOL);
    let restored = reader.body().unwrap(); reader.finish().unwrap();
    let roundtrip=write_body(&restored); assert_eq!(wire.len(),roundtrip.len());
    for (&a,&b) in wire.iter().zip(&roundtrip) {assert!((a-b).abs()<=1e-12*a.abs().max(1.));}
    let mesh = tessellate::tessellate(&restored,chord).unwrap();
    assert_eq!(wire,write_body(body));
    assert!(!mesh.triangles.is_empty());
    assert!(mesh.positions.iter().chain(&mesh.normals).all(|p| p.is_finite()));
    let bits = |x:f64| if x==0. {0} else {x.to_bits()};
    let key = |p:Vec3| (bits(p.x),bits(p.y),bits(p.z));
    let mut boundary = HashMap::<_,i32>::new();
    let mut volume = 0.;
    let origin = mesh.positions[0];
    for (tri,&face) in mesh.triangles.iter().zip(&mesh.faces) {
        let [a,b,c] = tri.map(|i|mesh.positions[i]);
        for (p,q) in [(a,b),(b,c),(c,a)] {
            if key(p)==key(q) {continue;}
            *boundary.entry((key(p),key(q))).or_default()+=1;
            *boundary.entry((key(q),key(p))).or_default()-=1;
        }
        let cross = b.sub(a).cross(c.sub(a));
        if cross.length()>1e-14 {
            let normal = tri.iter().fold(Vec3::new(0.,0.,0.), |s,&i|s.add(mesh.normals[i]));
            assert!(cross.dot(normal)>0.,"face {face}: inward triangle");
        }
        volume += a.sub(origin).dot(b.sub(origin).cross(c.sub(origin)))/6.;
        if let Surface::Cone {p,a:axis,r,k,..} = body.faces[face].surface {
            // The true distance to the conical surface bounds chord error at interior edge/face samples.
            for q in [a.add(b).scale(0.5),b.add(c).scale(0.5),c.add(a).scale(0.5),a.add(b).add(c).scale(1./3.)] {
                let w=q.sub(p); let h=w.dot(axis);
                let distance=(w.sub(axis.scale(h)).length()-(r+k*h)).abs()/(1.+k*k).sqrt();
                assert!(distance<=chord+1e-9,"face {face}: cone chord error {distance} exceeds {chord}");
            }
        }
    }
    assert!(boundary.values().all(|&n|n==0),"mesh is not closed across shared face edges");
    let measured=mass::measure(body,TOL).unwrap();
    assert!((volume-measured.volume).abs()<=chord*measured.area+1e-8);
    for (edge,_) in body.edges.iter().enumerate() {
        let owners:Vec<_>=body.faces.iter().enumerate().flat_map(|(f,face)|
            face.loops.iter().flatten().filter(move |&&(e,_)|e==edge).map(move |_|f)).collect();
        if owners.len()==2 && owners[0]==owners[1] {assert!(!mesh.display_edges.contains(&edge));}
    }
}

#[test]
fn tapered_solids_and_hollow_sections_have_analytic_mass_and_closed_outward_meshes() {
    for (rb,rt,h,wall) in [(5.,2.,30.,0.),(5.,2.,30.,0.25),(2.,5.,30.,0.1),(20.,0.5,2.,0.1),(2.,1.9,120.,0.01)] {
        let mut pts=vec![(rb,0.),(rt,h)];
        if wall==0. {pts.extend([(0.,h),(0.,0.)]);} else {pts.extend([(rt-wall,h),(rb-wall,0.)]);}
        let expected=frustum(rb,rt,h)-if wall==0. {0.} else {frustum(rb-wall,rt-wall,h)};
        let slant:f64=((rb-rt)*(rb-rt)+h*h).sqrt();
        let expected_area=PI*(rb+rt)*slant+PI*(rb*rb+rt*rt)
            +if wall==0. {0.} else {PI*(rb+rt-2.*wall)*slant-PI*((rb-wall).powi(2)+(rt-wall).powi(2))};
        for tilted in [false,true] {
            let mut s=section(pts.clone());
            if tilted {s.p=Vec3::new(17.,-23.,91.);s.a=Vec3::new(1.,2.,2.).unit().unwrap();s.x=Vec3::new(2.,-1.,0.).unit().unwrap();s.points.reverse();}
            let body=revolve(&s,TOL).unwrap();
            let m=mass::measure(&body,TOL).unwrap();
            assert!(m.exact);close(m.volume,expected);close(m.area,expected_area);
            checked_mesh(&body,0.01);
        }
    }
}

#[test]
fn stepped_flange_and_double_taper_are_general_sections_not_a_pole_special_case() {
    let variants = [
        (vec![(2.,0.),(2.,2.),(4.,2.),(4.,3.),(2.,3.),(2.,6.),(0.,6.),(0.,0.)],PI*(4.*5.+16.)),
        (vec![(5.,0.),(5.,3.),(3.,7.),(2.,9.),(2.,12.),(0.,12.),(0.,0.)],
         PI*25.*3.+frustum(5.,3.,4.)+frustum(3.,2.,2.)+PI*4.*3.)
    ];
    for (pts,expected) in variants {
        let body=revolve(&section(pts),TOL).unwrap();
        close(mass::measure(&body,TOL).unwrap().volume,expected);
        checked_mesh(&body,0.01);
    }
}

#[test]
fn direct_constructor_rejects_nonfinite_or_invalid_frames() {
    let valid=section(vec![(2.,0.),(1.,3.),(0.,3.),(0.,0.)]);
    for (axis,seam) in [(Vec3::new(f64::NAN,1.,0.),valid.x),(valid.a,Vec3::new(f64::NAN,0.,0.)),
                       (Vec3::new(0.,0.,0.),valid.x),(valid.a,valid.a)] {
        let mut s=valid.clone();s.a=axis;s.x=seam;
        assert!(revolve(&s,TOL).is_err(),"invalid frame admitted: {axis:?} {seam:?}");
    }
}

#[test]
fn direct_constructor_rejects_overflow() {
    let huge=section(vec![(1e200,0.),(1e200,1e200),(0.,1e200),(0.,0.)]);
    assert!(revolve(&huge,TOL).is_err(),"overflowing section admitted");
}

#[test]
fn wire_revolve_checks_lengths_and_preserves_exact_body() {
    let data=vec![0.,0.,0.,0.,1.,0.,1.,0.,0.,4.,2.,0.,1.,3.,0.,3.,0.,0.];
    let call=|data|protocol::handle(&protocol::Request {op:protocol::OP_REVOLVE,units:1,tolerance:TOL,data});
    let reply=call(data.clone()).unwrap();close(*reply.last().unwrap(),7.*PI);
    let mut reader=Reader::with_tolerance(&reply[..reply.len()-1],TOL);
    let body=reader.body().unwrap();reader.finish().unwrap();
    close(mass::measure(&body,TOL).unwrap().volume,7.*PI);
    assert!(call(data[..data.len()-1].to_vec()).is_err());
    let mut extra=data.clone();extra.push(90.);assert!(call(extra).is_err());
    for at in [0,3,6,10,11] {let mut bad=data.clone();bad[at]=f64::NAN;assert!(call(bad).is_err());}
}

#[test]
fn crossed_axis_apex_and_self_intersecting_profiles_fail_explicitly() {
    for pts in [vec![(-1.,0.),(2.,0.),(2.,3.),(-1.,3.)],vec![(0.,0.),(2.,0.),(0.,3.)],
        vec![(1.,0.),(3.,3.),(3.,0.),(1.,3.)],vec![(1.,0.),(1.,0.),(2.,2.)]] {
        assert!(revolve(&section(pts),TOL).is_err());
    }
}

#[test]
fn cylindrical_revolve_can_use_existing_through_cut_domain() {
    let body=revolve(&section(vec![(2.,0.),(2.,6.),(0.,6.),(0.,0.)]),TOL).unwrap();
    let tool=gc_kernel::extrude::extrude(&Profile {
        n:Vec3::new(0.,1.,0.),low:-1.,high:7.,loops:vec![vec![gc_kernel::geom::Curve::Arc {
            c:Vec3::new(0.,0.,0.),n:Vec3::new(0.,1.,0.),x:Vec3::new(1.,0.,0.),
            r:0.5,t0:0.,t1:2.*PI,
        }]],
    },TOL).unwrap();
    let cut=gc_kernel::boolean::cut(&body,&tool,TOL).unwrap();
    close(mass::measure(&cut,TOL).unwrap().volume,PI*(4.-0.25)*6.);
    checked_mesh(&cut,0.01);
}
