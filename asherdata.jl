using(DataFrames)
using(CSV)
using(RCall)
eigencluster = function(sp,evecs,n)
    carray = Array{Array}();
    carray[1] = sortperm(evecs[:,2]);
    println("Ordered List ","=",sp[carray[1]])
    for i=2:n
        eigenvec = evecs[:,i];
        rarray = Array{Array}(0);
        for j = 1:length(carray)
            set21 = find(x->x>0,eigenvec);
            set22 = find(x->x<0,eigenvec);
    
            set21 = intersect(carray[j],set21);
            set22 = intersect(carray[j],set22);
            if length(set21) > 0
                push!(rarray,set21);
            end
            if length(set22) > 0
                push!(rarray,set22)
            end
        end
        carray = copy(rarray);
    end
    for i=1:length(carray)
        println()
        println("Cluster ",i,"=",sp[carray[i]])
    end
    # println(showall(carray))
end




# [morphology matrix used in Asher RJ 2007 BMC Evolutionary Biology]


# DIMENSIONS  NTAX=53 NCHAR=196;
# FORMAT SYMBOLS= " 0 1 2 3 4 5" MISSING=? GAP=- ;
# CHARSTATELABELS  

# 1  stap_obturator_fmn_G96 / 0_large '1_small-abs',
# 2  'petrosal-squam joint' / 0_flat '1_ball(pet)-socket(sq)' 2_minimal_cntct 3_tubular,
# 3  squamosal_EAM / 0_no_approx_vent_EAM 1_squam_vent_EAM,
# 4  'fen rotund pres N89:73' / 0_present 1_absent,
# 5  'Fen rotund orient N89:42' / 0_posterolat 1_posterior,
# 6  'JugForamen W01:149' / 0_similar_FenCochleae 1_larger_FenCochleae,
# 7  'crista interfenestr W01:133' / 0_abs 1_pres,
# 8  width_basal_cochlea_S98 / '0_<20%_basicran_width' '1_>20%_basicran_width',
# 9  anterior_promontory / 0_ant_smooth 1_projects,
# 10  'Petrosal canals N89:81' / 0_absent 1_present,
# 11  'carotid sulcus W01:146' / 0_pres 1_abs,
# 12  'stapedial sulcus W01:147' / 0_pres 1_abs,
# 13  'pirif fenest M88:9,A99:2' / 0_absent 1_present,
# 14  petrosal_pachyostosis / 0_light 1_dense,
# 15  'tegtymp pneumat N89:59' / 0_laminar 1_air_cells,
# 16  mastoid_pneumat_CM86 / 0_laminar 1_air_cells,
# 17  'tegtymp process S94:6' / 0_flat 1_projects,
# 18  'epitympRec M88:10,A99:6' / 0_contin_midear 1_separate,
# 19  'basisphBulN86:49,M88:3,A99:1' / 0_reduced 1_present,
# 20  'Alisphenoid bulla N88:10' / 0_reduced 1_present,
# 21  basioccipital_bulla / 0_absent 1_shields_bulla_medially,
# 22  'CTPP N86:55,M88:2,A99:3' / 0_reduced 1_shields_postbulla,
# 23  'RTPP M88:1,A99:4' / 0_reduced 1_present,
# 24  'Ectotymp R98:196 A99:8' / 0_phaneric 1_occluded,
# 25  'ectotymShapeA99:8 W93:18' / 0_ringlike 1_expanded,
# 26  'ectotym-Medge T94,Lu98' / 0_abs 1_pres,
# 27  'extAud meatus -shape' / 0_ovoid 1_ventral_extension,
# 28  'tympanic annulus LG99:38' / 0_supports_membrane 1_conical_ligament,
# 29  'entotympanic N86:47,M88:7,A00:30' / '0_small-abs' 1_present,
# 30  'mastoidExp M94:12 N89:60' / 0_exposed 1_occluded,
# 31  'Pmast exp A99:17' / 0_posterior 1_lateral,
# 32  'forOvale N89:* G96' / 0_within_alisph 1_post_alisph 2_squamosal,
# 33  'alisph canal-PRS' / 0_pres 1_abs,
# 34  'alisph canal-REL' / 0_lat_V 1_confl_V,
# 35  'max-squam connection' / 0_unconnected 1_cnx_ventSphenOrb,
# 36  'ectopteryg N86:36,T92:23,A99:14' / 0_single 1_dual,
# 37  'basisphen pit A99:11' / 0_absent 1_present,
# 38  'postlacerateFm N86:67 T92:31' / 0_narrow 1_AP_elongate,
# 39  'hypoglossal foramen N89:53' / 0_abs 1_single 2_multiple,
# 40  ventral_basisphenoid / 0_exposed 1_covered_by_vomer 2_covered_by_pterygoids,
# 41  petmastoid_exp_braincase / 0_exposed 1_excluded,
# 42  'subarc fossa M94:15 N89:43' / 0_present 1_flat 2_sc_can_reduced,
# 43  'dors subarc N89:43' / 0_closed 1_dors_open,
# 44  Tent_cerebelli_W93_Sh98 / 0_absent 1_divides_brcase,
# 45  'sinus canal N86:38 W01:103' / 0_absent 1_present,
# 46  'dorsum sellae N86:37' / 0_weak 1_prominent,
# 47  'Crista galli N86,A99,Sh98' / 0_flat 1_dors_projects 2_cribriform_plate_abs,
# 48  'Fm.rotundum N86:33 M94:7' / 0_confluent_SF 1_distinct_SF,
# 49  alisphenoid_spine / 0_flat 1_anteriorly_projects,
# 50  'OptFm N89:25,A99:18,M94:1' / 0_absent 1_distinct,
# 51  'opticFm size A99:18 M94:1' / 0_caliber_as_V 1_caliber_smaller_V,
# 52  'suboptic fm A00:36,N89:13' / '0_absent-hidden' 1_present,
# 53  interorbital_fenestra / 0_absent 1_connects_orbits,
# 54  'ethmoid fmn T92:20' / 0_present 1_absent,
# 55  'TubMax N89:77,F93:67' / 0_absent 1_dorsal_ridge_in_orbit,
# 56  'Glenoid pos N89:51' / 0_even_petrosal 1_dors_petrosal 2_glenoid_undefined,
# 57  'glenoid shape N86:40' / 0_transverse 1_AP_elongate 2_dual,
# 58  postglenoid_foramen / 0_glenoid_imperforate 1_pgf_at_glenoid,
# 59  Postglenoid_process / 0_absent 1_present '2_pres-AP',
# 60  'ento-postglen processes' / 0_separate 1_connected,
# 61  'Entoglenoid proc A99:22' / '0_no-lat_jaw_support' 1_post_support,
# 62  'Zygo arch A99:26' / 0_complete 1_incomplete '2_squam-frontal',
# 63  'Jugal in glenoid N89:61' / 0_reaches_glenoid 1_reduced_fwd_glenoid,
# 64  'Anterior Jugal C95:15' / 0_antedge_orbit 1_lateral_orbit,
# 65  'Zygom proc squam N89:67' / 0_robust_dorslat 1_narrow_laterally,
# 66  'Postorb bar N89:80' / 0_reduced 1_surrounds_90%_orbit,
# 67  'Postorb septum R98:206' / 0_absent '1_jugal-alisph_cntct',
# 68  Interorb_distance / 0_separate 1_convergent,
# 69  facial_LacProc / '0_<=_orb_proc' 1_large_facial_proc,
# 70  'OrbLacrProc N89:78,F93:76' / 0_flat 1_proj_into_orbit,
# 71  Lacrimal_foramen / 0_present 1_not_expsd_antOrbit,
# 72  'LacForamen open N86:23' / 0_posterior 1_lateral,
# 73  'frontal in orbit N89:26' / 0_reduced_cntct_palatine 1_broad_cntct_palatine,
# 74  Jugal_in_MedOrbit / 0_jugal_restricted '1_jugal-parietal_cntct',
# 75  'Max medOrbit N86:14 A99:25' / 0_reduced 1_extensive,
# 76  'pmax-frontal contact N89:47' / 0_absent 1_present 2_pres_over_max,
# 77  anterior_nasals / 0_fused_pmax 1_projecting,
# 78  posterior_nasals / 0_tapered 1_broad 2_horns,
# 79  metopic_suture / 0_unfused 1_fused,
# 80  fenestrated_max / 0_solid_medZygoma 1_fenestr_medZygoma 2_trabeculated_medZyg,
# 81  InfOrb_canal_length / 0_long 1_short,
# 82  InfOrbCanal_antForamn / 0_small 1_similar_ext_nares,
# 83  'orbitPositn T92:17,N89:65' / 0_ant_alv_max 1_dors_alv_max 2_post_alv_max,
# 84  'externNares O99:39' / 0_ant_rostrum 1_dorsal_cheek_teeth 2_post_rostrum,
# 85  'postpalate spine N86:18' / 0_flat 1_prominent,
# 86  posterior_hard_palate / 0_post_alvmax 1_even_alvmax 2_ant_alvmax,
# 87  palatal_notch / 0_closed_medial_M3 1_open_medial_M3,
# 88  'postventral max A99:15' / 0_alvmax '1_postAlvMax-vent_pteryg',
# 89  'palatal pmax exp N89:50' / '0_abs-alveoli_only' 1_small 2_half_palate,
# 90  Incisive_fmn_presence / 0_present 1_absent,
# 91  'Incisive fmna S95:2 T92:18' / 0_dual_open_palate 1_single_opening_palate,
# 92  Incisive_Fmna_shape / 0_oval 1_elongate,
# 93  palatine_fenest / 0_solid 1_post_fenestr 2_ant_fenestr,
# 94  dental_reduction / 0_teeth_present 1_atrophied 2_baleen,
# 95  incisor_enamel / 0_covers_tooth 1_antBand '2_abs-tip_only',
# 96  'num incisors A99:38' / 0_none 1_one 2_two 3_three 4_four 5_five,
# 97  postPremax_diastema / '0_small-abs' 1_large_gap,
# 98  'incisor growth N89:54' / '0_crown-root_junct' 1_hypsodont,
# 99  'tooth replacement N89:11' / 0_antemolars 1_last_pm_only 2_uncalcified 3_continuous,
# 100  canine_growth_G01 / 0_determinate 1_hypsodont,
# 101  'cheektooth enamel N89:29' / 0_present 1_atrophied,
# 102  cheektooth_orientation / 0_parallel '1_V-shaped',
# 103  cheektooth_growth / '0_crown-root_suture' 1_hypsodont,
# 104  upper_P3 / 0_pres 1_abs,
# 105  'P3 crown A00:7' / 0_mediolat_compressed '1_prom_lingual_cusp(s)',
# 106  'P4 crown A00:71' / 0_mediolat_compressed '1_prom_lingual_cusp(s)',
# 107  'metacone A99:41' / 0_present 1_absent,
# 108  'protocone A99:42' / 0_present 1_absent,
# 109  molar_ectoloph / 0_unconnected 1_pi_shaped,
# 110  'stylarshelf N89:20,A00:79' / 0_broad 1_narrow,
# 111  'molar shearing W93:31' / 0_multilocus 1_carnassial,
# 112  'mandible M94:31' / 0_robust 1_threadlike,
# 113  condyle_height / 0_at_toothrow 1_above_toothrow,
# 114  coronoid_shape / 0_large 1_reduced,
# 115  coronoid_position / 0_superior_condyle 1_inferior_condyle,
# 116  jaw_angle / 0_narrow 1_medial_inflection,
# 117  'symphysis-fusion' / 0_unfused 1_fused,
# 118  'symphysis-length' / 0_AP_narrow '1_>30%_dentary_len',
# 119  'symphysis-vent lip' / 0_absent 1_anteroventral_lip,
# 120  'lower m3 presence A00:89' / 0_present 1_absent,
# 121  'talonid basin A99:44' / 0_present 1_reduced,
# 122  'cervCentra len O99:88' / '0_>_2x_thoracic' 1_similar_thoracic 2_AP_compressed,
# 123  C2_odontoid / 0_present 1_reduced,
# 124  C2_odontoid_shape / 0_broad '1_narrow-pointed',
# 125  'C3-7 spinous processes' / 0_present 1_reduced,
# 126  'ventrum of sternum ?S98:91' / 0_flat 1_keeled,
# 127  sternal_ribs / 0_unossified 1_ossified,
# 128  'Clavicle  A00:100' / 0_elongate 1_stout 2_absent,
# 129  'dorsum Proxribs N89:37' / 0_round 1_flattened,
# 130  'intervert artics M94:18' / 0_nomarthrous 1_xenarthrous,
# 131  'T5-T9 spinous processes' / 0_present 1_reduced,
# 132  'scapular shape S95:8' / 0_CC_long 1_ML_long 2_ML_long_and_narrow,
# 133  scapular_coracoid / 0_reduced 1_ventral_glenoid 2_anterior_glenoid,
# 134  'radius-humerus ratio' / 0_similar 1_radius_much_longer,
# 135  'HumeralHead M94:20' / 0_below_tuberos 1_even_tuberos 2_superior_tuberos,
# 136  'hum supinRidge S95:10' / 0_present 1_reduced,
# 137  'entepicond fm A99:57' / 0_present 1_absent,
# 138  'olecranon fossa A99:59' / 0_solid 1_fenestrated,
# 139  'med epicond A99:58' / 0_reduced 1_present 2_elongate,
# 140  dist_humerus_G01 / 0_cylindrical '1_V-shaped_troch',
# 141  proximal_radius / 0_capitulum '1_troch+capit',
# 142  prox_ulna / 0_single_surface 1_notched_for_rad,
# 143  'ulnar olec N89:41 S95:38' / '0_olec-humeral_artic' 1_reduced_olecranon,
# 144  'distUlnS95:38 A00:111 M94:21' / 0_no_carpal_artic 1_carpal_artic,
# 145  'ulna-radius shaft' / 0_separate 1_fused,
# 146  'distal radius-ulna' / 0_separate 1_fused,
# 147  metacarpal_length / 0_much_shorter_forearm '1_>=50%_forearm',
# 148  'Carpal arrangement N89:62' / 0_alternating 1_no_lunarUncif_art,
# 149  'scaph-lunateFusion W93:55' / 0_separate 1_fused,
# 150  '3rdPhalanx 2ndDigit S98:151' / 0_ossified 1_unossified,
# 151  term_phalanx_shape / 0_ML_narrow 1_ML_broad,
# 152  term_phalanx_number / '0_<=_three' '1_>_three',
# 153  'sacral verts A99:63' / 0_4_or_fewer '1_>4',
# 154  'pubicSymph M94:22,A99:62' / 0_broad 1_reduced,
# 155  ilium / 0_narrow 1_dorsal_process,
# 156  'IschNotch M94:23 (N89:69)' / 0_open 1_closed,
# 157  'sacral-ilium artic G98' / 0_present 1_absent,
# 158  ant_process_pubis / 0_absent 1_present,
# 159  'epipubic bones N89:19' / 0_present 1_absent,
# 160  'caudals distsacrum W93:54' / 0_elongate 1_reduced,
# 161  'fem fovea capitis M94:27' / 0_central 1_marginal 2_indistinct,
# 162  femoral_neck / '0_pres-head_med_shaft' '1_abs-head_above_shaft',
# 163  '3rdFemTrochA99:61 O99:101' / '0_pres-proximal' '1_pres-midshaft' 2_reduced,
# 164  med_troch_ridge / 0_even_lat 1_enlarged,
# 165  'distal tibia-fib A99:54' / 0_unfused 1_fused,
# 166  'fibular ossification S98:170' / 0_ossified 1_incomplete,
# 167  'fibfacet calc T92:46 M01:79' / '0_fib-calc_artic' '1_facet_small-abs',
# 168  'peroneal proc A99:71' / 0_blunt 1_large,
# 169  dorsum_astrag / 0_flat 1_trochlea,
# 170  'Sustentac facet N89:36' / 0_separate 1_cntct_dist_astrag,
# 171  'astrag cot fossa M94:28' / 0_abs 1_pres,
# 172  'astrag postmed proc M94:29' / 0_small_or_abs 1_prominent,
# 173  'astrag navicFacet O99:104' / 0_convex 1_saddle 2_trochlea 3_concave,
# 174  'calc facet of astrag O99:107' / 0_plantar 1_lateral,
# 175  'astrag-cuboid cntct G01' / 0_reduced_contact 1_contact,
# 176  prox_cuboid / 0_flat 1_stepped,
# 177  'calcar S94:19' / 0_absent 1_uropatag_spurGastr 2_present,
# 178  pedal_digital_ray_1 / 0_weight_bearing 1_reduced,
# 179  pedal_digital_ray_2 / 0_weight_bearing 1_reduced,
# 180  pedal_digital_ray_3 / 0_weight_bearing 1_reduced,
# 181  pedal_digital_ray_4 / 0_weight_bearing 1_reduced,
# 182  pedal_digital_ray_5 / 0_weight_bearing 1_reduced,
# 183  1st_metatarsal / 0_as_other_digits 1_opposable,
# 184  metatarsal_length / 0_shorter_tibia 1_approach_50%_tibia,
# 185  'hallux proxPhal S94:21' / 0_as_other_digits 1_elongate,
# 186  forelimb_patagium / 0_absent 1_present,
# 187  'cloaca (ShMcK98:215)' / '0_anal-urogen_separate' 1_present,
# 188  trophoblast / 0_abs 1_pres,
# 189  'placenta (Carter, 2001)' / 0_choriovitelline 1_chorioallantoic,
# 190  'placenta2 (Carter, 2001)' / 0_epitheliochorial 1_endotheliochorial 2_haemochorial,
# 191  'yolk sac (Carter, 2001)' / 0_permanent 1_temporary 2_rudimentary,
# 192  gestation_time / 0_short 1_prolonged,
# 193  'ureter-muellerian ducts' / 0_?? 1_lateral,
# 194  hindgut / 0_differentiated 1_simplified,
# 195  'testes (werd-nils 1999)' / 0_abdominal '1_descend-ascr' '2_descend-scr',
# 196  eustacian_sac / 0_abs 1_pres



sp = [
"Didelphis",
"diprotodont",
"Amblysomus",
"balaenopterid",
"Caniform",
"Cavia",
"Ceratotherium",
"Chaetophractus",
"Cynocephalus",
"delphinid",
"Echinops",
"Elephantulus",
"Equus",
"Erinaceus",
"Felis",
"Hippopotamus",
"Homo",
"Hystrix",
"Lama",
"leporid",
"Loxodonta",
"Macroscelides",
"Manis",
"Mus",
"myrmecophagid",
"Ochotona",
"Orycteropus",
"phyllostomid",
"Procavia",
"Pteropus",
"Rattus",
"Rousettus",
"ruminant",
"Solenodon",
"Sorex",
"strepsirhine",
"Sus",
"Talpa",
"Tapirus",
"Trichechus",
"Tupaia",
"Anagale",
"Arsinoitherium",
"Centetodon",
"Hyopsodus",
"Leptictis",
"Meniscotherium",
"Moeritherium",
"Palaeoparadoxia",
"Phenacodus",
"Plesiorycteropus",
"Ukhaatherium",
"Zalambdalestes"
];
nsp = length(sp);

charstate = Array{String}(nsp);

charstate[1] =           "(01)000010000110000000100000-00000(01)1-00002000000001-0-(01)00000010000000-0000010000100001000001001100500100000000000000001000001010000000100100010000100000000000000002020001(01)010100000000001000100--00020";
charstate[2] =         "000001000011000001010000100001-01-0100200000000100-010(01)00010000000-0000110000100001000001001(01)013101-0001-10001001001000001010000000100(01)00010000100001000000000000020001000100001011100-1-0100--00020";
charstate[3] =          "0000011001000011-1110110100000011-00001000001000-0-00001010-00--10-0?001?-100?101010(01)000100000030000000011100000011100011101110000010020002000010000010001000011201010010000000000000000000112011100";
charstate[4] =        "13001100001111--00111000110100021-01010112-00020-0-001000010021010-00-1-1002-0000022000001--02-0---------------00100000--21-00?2000000211001-01100000001----1-10-------------------------00110111000";
charstate[5] =             "00000100001100000000000010001000(01)(01)00001000000101010(01)00000010001010-0000010000000001001(01)01001000300000000000001100000000001010002000100011100000100001000000000100020001010000000010000-0-00111011020";
charstate[6] =                "0000010010110000000000001010001(01)0011012000000000-10110111020000110-00100001100011110020020010011110-0111-10001001110000001010002000100011100100100011000000000110010001011000000010001-0-00112011010";
charstate[7] =        "000001000011100010000100000001-10100011001000010-1100?100010001000-001001000111000100200001000-00-0-00101100110010001000010(01)0102100100001000100100000010001000111011001010001010010001-0-00110111010";
charstate[8] =       "000001000011000111000000100011-01-00001001001010-10000010011000010-000000010010010100000100000-10-01101-------001000000--101101011011010001010010000000011010010101010001000000001000100000112011010";
charstate[9] =         "(01)010111000(01)1000101000000100010001-00011000000010-1000000011000011100000010000110001012001001000(01)000000000100000001000000010(01)001010011120011010100-00100001000010000000101100000000000000010112?110?0";
charstate[10] =            "12001100001111--00101000110101-01-01011112-00020-0-001000010001010000-1-1012-1000022000011--000200200000001101000100000-121-1012000020111001-01100000011----1-10-------------------------00110111100";
charstate[11] =            "0000011000000000001100000-0000101-0010(12)000001000-1100000000-11--10-000011-1000001010000110000002000000001111000000000001110110000001100000100001000000000100001101(02)000111100000000000000001112011110";
charstate[12] =         "00000110011000000001011010001010(01)001001000001000-1011011100-001010-0000010000000101011101001100300000000110001001010000(01)010110000001101101001000--00000000000010000010001011000001000001000112011000";
charstate[13] =                "000001000011100110000100100000010100011001000011010000000010001001000(01)00100011(01)00010020010010003100000101100110010001000000010020001000010001000--100-10101000111011-1-010001000011011-1-00110111021";
charstate[14] =            "000001100000(01)001001000100-000000(01)001101000001000-11000000110001110-000000010000000101000100010030000000001000100000000000101100000011011110010010000000001000010201010001000000001000000000112011110";
charstate[15] =                "000001000011000000000000100010001-0000100001010101000000001000101101000010000000101000101000000300000000000001100000000101010002000100(01)10000000100001000000000100020001010000000010000-0-00111011020";
charstate[16] =         "0100010000111001?0000000100001-11-00011001000010-10000100010001011001010000001110010000010000002110100000100010000001(01)0001010102000100001000100111000010101000112021000010002111010000-0-00110111020";
charstate[17] =                 "0000011001110011-10101101000001(01)1-010110010001110100000000100011111100000000001000000(01)00101000020000000011000100100010000101000010001021101000010000000000100011002100101100000000000000000112111020";
charstate[18] =              "000001001011000001000000100001-10001001000000010-10000011020000110-00100001101001110020020000011110-0011-100010011100100010100000001100(01)111010010000100000000011002000101100000001000000000112011010";
charstate[19] =                "000001?000??100000000000100001-(01)1-01012000000010-1000010011000101100010000000100001002101000000110000000010001001000100000000002000100011000100111100-10(01)01000110020010010002111011001-1-00110111020";
charstate[20] =              "000001000001000100000000100000110001002000000110-1001011000-000?1100--00001101021010(01)20020012012110-0010110001001110000001011102000110011100100100000000000000110000100010000000010000-1-00112011020";
charstate[21] =           "0011---000?10111-0000000100011-1010000000100001111000010000-000100-0011-00?1-1101001020010110021113-000---000100101010100201010200010010100001010001101000100010202000000001000000000000000111111000";
charstate[22] =       "000010101100000101110110100010(01)01-01001000001000-1011011000-000010-0000010000000101011001001(12)003000001001100010010100001010110000001101101001000--00000000000010000010001011001001000001000111011000";
charstate[23] =               "100001000001000001000000100011-01-00001001010(01)10-1100000-00-001110-0--(01)01000011010200000100001-0--2------------1-110000--1010002000100200010100100001000010(01)0010202000011000301001000100000110011110";
charstate[24] =                 "00000100001010000000000010000010000100(12)00000(01)000-1001011100-001110-0--1-000100011110000020011011110-0001-100010011110000010110000011102011101001000010000100001000(12)010111000000000000000000112011020";
charstate[25] =       "000001000001000000011110100011-01-0000120000001(01)01000000100-011110-00001100001000020-00001--01-0--2------------1-110000--101001211012020001000010000000011010010012000101000300000000000000112?11010";
charstate[26] =             "0000110001110011-000011010000011000100(01)000000000-1001011000-000110-00000000100011010021020(01)(01)2012110-0010110001001110000001011100000110111100100100000000010000110000100010000000010000-0-01112011020";
charstate[27] =         "0000010000000000010000000-0000001-00001001011110-10000100011001010-010001000010000100000100000-01-0-101-------001000000--1010000000120100(01)1000010000001010100110101000111011000001000100000111211010";
charstate[28] =        "000011010011100000000000100000000100(01)11000101101111000000010101?(01)0-0?0(01)0?0?0??10001000(01)010(01)00002000000001100010000001000010111(01)01010110110101010--10110001000111010001100100010020000000110112011020";
charstate[29] =            "000001000011100000010000100011-00101011000000011(01)100(01)010001010000100010010000100101012001000000101000000110011001000100001010002000100011100100100010010100000111000001010100000010001-0-00112111001";
charstate[30] =             "000011000010(01)0001000000010000000(01)101011000000010-1000000001000111(01)000001100000101010000000100002000000001100010000001000010111101010111110101010--10100001010111012001100100011010000000110112011010";
charstate[31] =               "000001100010000000000000100000100101001000000000-1001011100-001110-0000-000100111110000020010011110-0001-10001001100000001010(01)00000110(12)01(01)101001000010000(01)000010001010111000000000000000000112011020";
charstate[32] =            "000011000000100010000000100000000100011000001000-10000000010001110-0?001?00000101010000000100002000000000100010000001000010111001010111110101010--10100001010111(02)1200110?100011010000000110112011010";
charstate[33] =             "000001000010000100000000100000(01)01-00012001010110-100101000101010110010000000000000101210100100-01-0-0000110001001000000000000002000100011000100111100-00(01)010001(01)0021010011002111011001-1-00110111020";
charstate[34] =            "000001100000100010000110000000100100012000001000-1100000000-11--10-0?0000010001010100000100000030000000011100000000001001101001000010010011000010000000011000010000000101000001000000000000111011110";
charstate[35] =                "0000010000001000000001000-0000101-0000(01)000101000-0-00000200-11--10-0-001?-1000001010000010000003002000000100000000000000010110000012001000100001000010000100001000001011?000000001000001000111011110";
charstate[36] =         "0000010001010000010001110-0000001-0100100001011(01)0100000000100010110100011000000000100200100000020000000011000100000000000101000000011020001000010000000000000010000000101110000000000010000110111020";
charstate[37] =                  "000001000001000000000000100001-11-01011000000110-1000000001000000(01)00110100000110002000001001000300010000110001001010100001010102000100001000100100100010001000110020000010002111011001-1-00110111020";
charstate[38] =                "0000011001000010-0110110100000100100001000(01)01000-(01)100001000-00--10-0?001?-100?0010100100100010030000000001000000100000000101110100120000001000010000000011000011100010111000001000000000000110011110";
charstate[39] =              "000001000011101110000000100001-10100011001000011010000100011001000-0010(01)10001100101112001011000310000000110011001000100001000102000100011(01)00100100100010101000111010001010001010010001-0-00110111011";
charstate[40] =           "1101---0001111--000000000-0001-11-00011001000010-1100110001000000100--1-1001--1111010200101100-01-3-0000110001001000101-020110021001001110101001010100100-001010-------------------------00111111000";
charstate[41] =               "0000010001000000000001010-0010110001001000001000-100100000100001110001001000011000100100100110020000000001000000100000000101000000011020001000010001100000000010000000101101000000000000000111011020";
charstate[42] =              "00001100010000000100011010001010000100200000?1?0-10010000110000?10-00101?0000100001011101???100????00000110001001000000001010???0001?0200?10?00??????000????0??000?000??10?0?0?0?00000?0????????????";
charstate[43] =       "?011---0001110111000000???0??1-10000010001?????0-10000100000001000-0001-101002?0101002001010000300000000110001001000110002011???100000101000110100010???0010001?112000000000001000000000????????????";
charstate[44] =           "?00001?000011000000001100000001001000?10000010?0-1100000000-10?110-0000000100000101001001000000300000000110000000000000001011??????11??00?1??00?????????????????????10??????????????????????????????";
charstate[45] =            "?000(01)100000100110000000???0??0000100011000??0??0-10001000010001000-00000000001000010020010000003000000001100010010001000010100?0000120210110100100000?000?0?00100010001110000000000000000???????????";
charstate[46] =            "000001000000(01)000000000110-0010000001001000001010-10000000010001010-001001000000010101(12)0010000002000000001100010000000000010111??000100200110100100000000000000100000101010000000000000000???????????";
charstate[47] =       "?000(01)1000000101?0000000???0??000000(01)011000?????1?1000?00001000100(01)000(01)0010?00000001002001000000300000000110001001000100001010002000120000110100100000010000000101010000010110000010001000???????????";
charstate[48] =         "??11---00011?0?00110000???0??1-000000?000100???101000000000-000100-0?01-00100110100001001010000311?0000011000100101000000?010?????0??0011010?101????????0?0?001?001100??0???????????????????????????";
charstate[49] =      "?01???????????????000?????0??1-0??000110?????????10?0?100010?0110100001-?0?0010010011000201100(02)3113000001100010010001100010100?200010001100010010000001000000011?020001010000110010000-0-???????????";
charstate[50] =           "?00001000010101?0000000???0??0000(01)0001(12)000??1??0-10000000010000010-001001000010010100(12)001000000300000000110001001000000001010??000010011011010010000001000000010101(01)001010000000010001000???????????";
charstate[51] =     "?000011000000000000000(01)???0??0001-0?001001000??101100??0000-101?10-0?100?0?0001???10?????????1??????-????????????????????????????0???0200110000100????0?11000110101010??0011000????????0????????????";
charstate[52] =         "?000001000101000000100001000?00?1-000020000?1????1??0?000111001010-00000?0?00100001000001000000500?00000010000000000000001??1??0000110100010000100000?000000000??02000110100000000000??0????????????";
charstate[53] =      "?000001000101000100000001000?0011-010020000010?1010?00100111001010-00000100001000010010?1000001(23)11?000000100010010000000010110?0000110010??01001000???0000000000?020101011000000010000010???????????";

#clean up the data
# (01) --> (
charmod = Array{String}(nsp);
for i = 1:nsp
    insertp = find(x->x=='(',charstate[i]);
    if length(insertp) > 0
        deletechar = Array{Int64}(length(insertp)*3)
        t = 1;
        for k=1:length(insertp)
            deletechar[t] = insertp[k];
            deletechar[t+1:t+2] = collect(insertp[k]+1:insertp[k]+2);
            t = t+3;
        end
        keep = deleteat!(collect(1:length(charstate[i])),deletechar);
        charmod[i] = charstate[i][keep];
    else
        charmod[i] = charstate[i];
    end
end

#Build character matrix
nmeas = length(charmod[1]);
CM = Array{Char}(nsp,nmeas);
for i = 1:nsp
    for j=1:nmeas
        CM[i,j] = charmod[i][j];
    end
end
#I think we need to treat ? and - the same
#Build similarity matrix
PC = Array{Float64}(nsp,nsp);
# measmeans = mean(pcdatatr[!isnothing(pcdatatr)],1);
#Build similarity matrix
for i = 0:(nsp^2 - 1)
    a = mod(i,nsp) + 1;
    b = Int64(floor(i/nsp)) + 1;
    
    if a == b
        PC[a,b] = 0.0;
        continue
    end
    # if a==1
    #     println(b)
    # end
    ct = 0;
    ctones = 0;
    for j = 1:nmeas
        if CM[a,j] != '?' && CM[b,j] != '?' 
            if CM[a,j] != '-' && CM[b,j] != '-'
                if CM[a,j] == '(' || CM[b,j] == '('
                    ct += 1;
                else
                    if CM[a,j] == CM[b,j]
                        # ct += (sqrt((pcdatatr[a,j]-pcdatatr[b,j])^2)); #/mean(pcdatatr[!ismissing.(pcdatatr[:,j]),j]);
                        ct += 1;
                    end
                end
                ctones += 1;
            end
        end
        # ctones += PA[a,j] + PA[b,j];
    end
    PC[a,b] = Float64(ct/ctones); #/Float64(ctones);
end


S = zeros(Float64,nsp,nsp);
# S = copy(-PC);
for i = 1:nsp

    val = zeros(Float64,10);
    lok = zeros(Int64,10);

    for j = 1:nsp
        if PC[i,j] > val[10]
            val[10] = PC[i,j];
            lok[10] = j;
        end
        for k=1:9
            if val[11-k] > val[10-k]
                v = val[11-k];
                val[11-k] = val[10-k];
                val[10-k] = v;
                l = lok[11-k];
                lok[11-k] = lok[10-k];
                lok[10-k] = l;
            end
        end
    end
    S[i,lok] = -PC[i,lok];
    S[lok,i] = -PC[lok,i];
end
rowsums = sum(S,2);
S[diagind(S)] = -rowsums;

ev = eigs(S; nev=10,which=:SR);
eval = ev[1];
evecs = ev[2];

ranked = sortperm(evecs[:,2]);
sp[sortperm(evecs[:,2])]

R"""
image($(PC[ranked,ranked]))
"""

R"""
plot($(evecs[:,2]),$(evecs[:,3]))
text($(evecs[:,2]),$(evecs[:,3]),$sp,cex=0.5)
"""

R"""
library(plot3D)
scatter3D($(evecs[:,2]),$(evecs[:,3]),$(evecs[:,4]))
text3D($(evecs[:,2]),$(evecs[:,3]),$(evecs[:,4]),labels=$sp,cex=0.5,add=T)
"""
    
