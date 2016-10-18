# m=1 reflection grating 

def get_beamline(ekev):
    import os
    import wpg
    from wpg import Beamline
    from wpg.optical_elements import Aperture, Drift, Lens, Empty, Use_PP, Mirror_plane, VLS_grating
    from wpg.wpg_uti_oe import show_transmission
    from wpg.srwlib import SRWLOptMirPl
    import numpy as np

    wpg_path = os.path.abspath(os.path.dirname(wpg.__file__))
    data_path = 'data_example_04'

    # S1 beamline layout
    # Geometry ###
    src_to_hom1 = 268.8 #274.  # Distance source to HOM 1 [m]
    src_to_hom2 = 271.7 #284.  # Distance source to HOM 2 [m]
    
    src_to_m3a = 287.4    # high energy 1-3 keV #Distance source to vertical focusing mirror of grating mono
    src_to_m3b = 288.2    # low energy 0.27 1.2 keV #Distance source to vertical focusing mirror of grating mono
    src_to_vls_gr = 288.8 #
    src_to_m5 = 326.8   # SCS distribution mirror
    src_to_pslit = 387.9  #400. # distance source to slit in focus 
    
    src_to_m3 = src_to_m3b #low energy mono

    theta_om = 9.e-3   # [rad] check!
    theta_m3 = 20.e-3  # [rad]
    alpha_gr1 = 16.239e-3 # [rad] for 0.8 keV
    vls_gr1_par = [50., 0.00101, 0, 0] #Polynomial coefficients for VLS Grating Groove Density
    vls_gr2_par = [150., 0.00101, 0, 0] #Polynomial coefficients for VLS Grating Groove Density
    # calculate grating geometry
    d = 1.e-3/vls_gr1_par[0]; #period
    wl = 12.39e-10/ekev
    m = 1 #grating order
    dtheta = np.arcsin(m*wl/(2*d*np.sin(theta_m3))) #for constant deflection angle geo
    alpha = theta_m3-dtheta; beta = theta_m3+dtheta
    print('alpha, beta: {:.1f}, {:.1f} mrad'.format(alpha*1e3,beta*1e3))
    print('{:.4g}'.format(np.cos(alpha)-np.cos(beta)))
    print('{:.4g}'.format(m*wl/d))    
    om_mirror_length = 0.8  # [m]
    om_clear_ap = om_mirror_length * theta_om
    
    m3_mirror_length = 0.58
    m3_clear_ap = m3_mirror_length * theta_m3
    gr_length = 0.5
    gr_width = 0.2
    gr_clear_ap = gr_length * theta_m3
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
    pm3_wavefront_distortion = Mirror_plane(orient='y', 
                                             theta=theta_m3, 
                                             length=m3_mirror_length, 
                                             range_xy=m3_clear_ap/aperture_x_to_y_ratio, 
                                             filename=os.path.join(
                                             mirrors_path, 'mj37_2.dat'), 
                                             scale=5.,
                                             bPlot=True)
    print('M3 WF distortion'); show_transmission(pm3_wavefront_distortion);
    bl0.append(pm3_wavefront_distortion,
               Use_PP(semi_analytical_treatment=0))
    q_m3 = src_to_pslit - src_to_m3
    bl0.append(Lens(1e23,1./(1./src_to_m3 + 1./q_m3)), 
            Use_PP(semi_analytical_treatment=1))
    m3_to_vls_gr = Drift(src_to_vls_gr - src_to_m3)
    bl0.append(m3_to_vls_gr, Use_PP(semi_analytical_treatment=1))
    vls_gr_0 = Aperture('r', 'a', m3_clear_ap, m3_clear_ap / aperture_x_to_y_ratio)
    #Grating Substrate (plane mirror, deflecting in vertical plane):
    gr_sub = SRWLOptMirPl(_size_tang=gr_length, _size_sag=gr_width, _ap_shape='r',
                    _nvx=0, _nvy=np.cos(alpha_gr1), _nvz=-np.sin(alpha_gr1), _tvx=0, _tvy=np.sin(alpha_gr1))    
    vls_gr_1 = VLS_grating(_mirSub=gr_sub, 
                           _m=1, 
                           _grDen=vls_gr1_par[0], 
                           _grDen1=vls_gr1_par[1], 
                           _grDen2=vls_gr1_par[2], 
                           _grDen3=vls_gr1_par[3])
    
    bl0.append(vls_gr_0, Use_PP(semi_analytical_treatment=0,
                            zoom=zoom, sampling=zoom / 1.))
    bl0.append(vls_gr_1, Use_PP(semi_analytical_treatment=0,
                            zoom=zoom, sampling=zoom / 1.))
    vls_gr_to_foc = Drift(src_to_pslit - src_to_vls_gr)
    bl0.append(vls_gr_to_foc, Use_PP(semi_analytical_treatment=1))
    bl0.append(Empty(),Use_PP(zoom_v=0.1, sampling_v=0.1/0.1))
    return bl0