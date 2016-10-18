# zero order (mirror) reflection

def get_beamline():
    import os
    import wpg
    from wpg import Beamline
    from wpg.optical_elements import Aperture, Drift, Lens, Empty, Use_PP, Mirror_plane, VLS_Grating
    from wpg.wpg_uti_oe import show_transmission

    wpg_path = os.path.abspath(os.path.dirname(wpg.__file__))
    data_path = 'data_example_04'

    # S1 beamline layout
    # Geometry ###
    src_to_hom1 = 274.  # Distance source to HOM 1 [m]
    src_to_hom2 = 284.  # Distance source to HOM 2 [m]
    src_to_m3 = 300.    # Distance source to vertical focusing mirror of grating mono
    src_to_vls_gr = 301.
    src_to_pslit = 400. # distance source to slit in focus 

    theta_om = 9.e-3  # [rad]
    theta_m3 = 9.e-3  # [rad]

    om_mirror_length = 0.8  # [m]
    om_clear_ap = om_mirror_length * theta_om
    
    m3_mirror_length = 0.5
    m3_clear_ap = m3_mirror_length * theta_m3
    # define the beamline:
    bl0 = Beamline()

    # Define HOM1.
    aperture_x_to_y_ratio = 1
    hom1 = Aperture(
        shape='r', ap_or_ob='a', Dx=om_clear_ap, Dy=om_clear_ap / aperture_x_to_y_ratio)
    bl0.append(
        hom1, Use_PP(semi_analytical_treatment=0, zoom=1,sampling = 1/0.6))

    # Define mirror profile
    # Apply distortions.
    mirrors_path = os.path.join(wpg_path, '..','samples', 'data_common')
    hom1_wavefront_distortion = Mirror_plane(orient='x', 
                                             theta=theta_om, 
                                             length=om_mirror_length, 
                                             range_xy=om_clear_ap/aperture_x_to_y_ratio, 
                                             filename=os.path.join(mirrors_path, 'mirror1.dat'), 
                                             scale=1, delim=' ',bPlot=True)
    print('HOM1 WF distortion'); show_transmission(hom1_wavefront_distortion);

    zoom = 1.2

    bl0.append(hom1_wavefront_distortion,
               Use_PP(semi_analytical_treatment=0, zoom=zoom, sampling=zoom/0.8))

    # Free space propagation from hom1 to hom2
    hom1_to_hom2_drift = Drift(src_to_hom2 - src_to_hom1)
    bl0.append(hom1_to_hom2_drift, Use_PP(semi_analytical_treatment=0))

    # Define HOM2 as aperture.
    zoom = 1.0
    hom2 = Aperture('r', 'a', om_clear_ap, om_clear_ap / aperture_x_to_y_ratio)
    bl0.append(hom2, Use_PP(semi_analytical_treatment=0,
                            zoom=zoom, sampling=zoom / 1))

    # drift to M3 mirror
    hom2_to_m3 = Drift(src_to_m3 - src_to_hom2)

    bl0.append(hom2_to_m3, Use_PP(semi_analytical_treatment=1))
    m3 = Aperture('r', 'a', m3_clear_ap, m3_clear_ap / aperture_x_to_y_ratio)
    bl0.append(m3, Use_PP(semi_analytical_treatment=0,
                            zoom=zoom, sampling=zoom / 1))
    # Define mirror profile
    # Apply distortions.
    mirrors_path = 'data_common'
    # Define mirror profile
    mm3_wavefront_distortion = Mirror_plane(orient='y', 
                                             theta=theta_m3, 
                                             length=m3_mirror_length, 
                                             range_xy=m3_clear_ap/aperture_x_to_y_ratio, 
                                             filename=os.path.join(
                                             mirrors_path, 'mj37_2.dat'), 
                                             scale=5.,
                                             bPlot=True)
    print('M3 WF distortion'); show_transmission(mm3_wavefront_distortion);
    bl0.append(mm3_wavefront_distortion,
               Use_PP(semi_analytical_treatment=0))
    q_m3 = src_to_pslit - src_to_m3
    bl0.append(Lens(1e23,1./(1./src_to_m3 + 1./q_m3)), 
            Use_PP(semi_analytical_treatment=1))
    vls_gr_0 = Aperture('r', 'a', m3_clear_ap, m3_clear_ap / aperture_x_to_y_ratio)
    bl0.append(vls_gr_0, Use_PP(semi_analytical_treatment=0,
                            zoom=zoom, sampling=zoom / 1.))
    vls_gr_to_foc = Drift(src_to_pslit - src_to_vls_gr)
    bl0.append(vls_gr_to_foc, Use_PP(semi_analytical_treatment=1))
    bl0.append(Empty(),Use_PP(zoom_v=0.1, sampling_v=0.1/0.1))
    return bl0