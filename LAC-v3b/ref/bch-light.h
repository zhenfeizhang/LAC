uint16_t a_pow_tab[256] = {1,2,4,8,16,32,64,128,29,58,116,232,205,135,19,38,76,152,45,90,180,117,234,201,143,3,6,12,24,48,96,192,157,39,78,156,37,74,148,53,106,212,181,119,238,193,159,35,70,140,5,10,20,40,80,160,93,186,105,210,185,111,222,161,95,190,97,194,153,47,94,188,101,202,137,15,30,60,120,240,253,231,211,187,107,214,177,127,254,225,223,163,91,182,113,226,217,175,67,134,17,34,68,136,13,26,52,104,208,189,103,206,129,31,62,124,248,237,199,147,59,118,236,197,151,51,102,204,133,23,46,92,184,109,218,169,79,158,33,66,132,21,42,84,168,77,154,41,82,164,85,170,73,146,57,114,228,213,183,115,230,209,191,99,198,145,63,126,252,229,215,179,123,246,241,255,227,219,171,75,150,49,98,196,149,55,110,220,165,87,174,65,130,25,50,100,200,141,7,14,28,56,112,224,221,167,83,166,81,162,89,178,121,242,249,239,195,155,43,86,172,69,138,9,18,36,72,144,61,122,244,245,247,243,251,235,203,139,11,22,44,88,176,125,250,233,207,131,27,54,108,216,173,71,142,1};

uint16_t a_log_tab[256] = {0,0,1,25,2,50,26,198,3,223,51,238,27,104,199,75,4,100,224,14,52,141,239,129,28,193,105,248,200,8,76,113,5,138,101,47,225,36,15,33,53,147,142,218,240,18,130,69,29,181,194,125,106,39,249,185,201,154,9,120,77,228,114,166,6,191,139,98,102,221,48,253,226,152,37,179,16,145,34,136,54,208,148,206,143,150,219,189,241,210,19,92,131,56,70,64,30,66,182,163,195,72,126,110,107,58,40,84,250,133,186,61,202,94,155,159,10,21,121,43,78,212,229,172,115,243,167,87,7,112,192,247,140,128,99,13,103,74,222,237,49,197,254,24,227,165,153,119,38,184,180,124,17,68,146,217,35,32,137,46,55,63,209,91,149,188,207,205,144,135,151,178,220,252,190,97,242,86,211,171,20,42,93,158,132,60,57,83,71,109,65,162,31,45,67,216,183,123,164,118,196,23,73,236,127,12,111,246,108,161,59,82,41,157,85,170,251,96,134,177,187,204,62,90,203,89,95,176,156,169,160,81,11,245,22,235,122,117,44,215,79,174,213,233,230,231,173,232,116,214,244,234,168,80,88,175};

uint32_t mod8_tab[1024] = {0,486539264,973078528,654311424,1946157056,1761607680,1308622848,1392508928,3892314112,4110417920,3523215360,3472883712,2617245696,2164260864,2785017856,3137339392,3439329280,
3489660928,4143972352,3925868544,3103784960,2751463424,2197815296,2650800128,620756992,939524096,520093696,33554432,1358954496,1275068416,1795162112,1979711488,2264924160,
2583691264,3170893824,2684354560,4076863488,3992977408,3372220416,3556769792,1862270976,1912602624,1426063360,1207959552,452984832,100663296,553648128,1006632960,1241513984,
1459617792,1879048192,1828716544,1040187392,587202560,67108864,419430400,2717908992,3204448256,2550136832,2231369728,3590324224,3405774848,3959422976,4043309056,318767104,
234881024,687865856,872415232,1728053248,2046820352,1560281088,1073741824,4211081216,3858759680,3238002688,3690987520,2399141888,2449473536,3036676096,2818572288,3724541952,
3271557120,3825205248,4177526784,2852126720,3070230528,2415919104,2365587456,905969664,721420288,201326592,285212672,1107296256,1593835520,2013265920,1694498816,2483027968,
2298478592,2919235584,3003121664,3758096384,4244635648,3657433088,3338665984,2080374784,1627389952,1174405120,1526726656,134217728,352321536,838860800,788529152,1493172224,
1140850688,1660944384,2113929216,754974720,805306368,385875968,167772160,2969567232,2885681152,2332033024,2516582400,3305111552,3623878656,4278190080,3791650816,637534208,
989855744,469762048,16777216,1375731712,1325400064,1744830464,1962934272,3456106496,3539992576,4093640704,3909091328,3120562176,2801795072,2147483648,2634022912,3942645760,
4127195136,3506438144,3422552064,2667577344,2181038080,2768240640,3087007744,50331648,503316480,956301312,603979776,1996488704,1778384896,1291845632,1342177280,2701131776,
3154116608,2600468480,2248146944,3573547008,3355443200,4009754624,4060086272,1224736768,1409286144,1929379840,1845493760,1023410176,536870912,117440512,436207616,1811939328,
1895825408,1442840576,1258291200,402653184,83886080,570425344,1056964608,2214592512,2566914048,3187671040,2734686208,4026531840,3976200192,3388997632,3607101440,889192448,
671088640,251658240,301989888,1090519040,1543503872,2063597568,1711276032,3707764736,3221225472,3875536896,4194304000,2835349504,3019898880,2466250752,2382364672,4160749568,
3841982464,3254779904,3741319168,2348810240,2432696320,3053453312,2868903936,268435456,218103808,704643072,922746880,1677721600,2030043136,1577058304,1124073472,2986344448,
2936012800,2281701376,2499805184,3321888768,3674210304,4227858432,3774873600,1509949440,1191182336,1610612736,2097152000,771751936,855638016,335544320,150994944,2130706432,
1644167168,1157627904,1476395008,184549376,369098752,822083584,738197504,2533359616,2315255808,2902458368,2952790016,3808428032,4261412864,3640655872,3288334336,0,
1275068416,2550136832,3556769792,754974720,1627389952,3036676096,4177526784,1509949440,369098752,3254779904,2382364672,1996488704,989855744,4009754624,2734686208,3019898880,
4160749568,738197504,1610612736,2566914048,3573547008,16777216,1291845632,3992977408,2717908992,1979711488,973078528,3271557120,2399141888,1526726656,385875968,1962934272,
956301312,3976200192,2701131776,1476395008,335544320,3221225472,2348810240,788529152,1660944384,3070230528,4211081216,33554432,1308622848,2583691264,3590324224,3238002688,
2365587456,1493172224,352321536,3959422976,2684354560,1946157056,939524096,2600468480,3607101440,50331648,1325400064,3053453312,4194304000,771751936,1644167168,3925868544,
2785017856,1912602624,1040187392,3338665984,2332033024,1593835520,318767104,2952790016,4227858432,671088640,1677721600,2634022912,3506438144,83886080,1224736768,1577058304,
301989888,3321888768,2315255808,1929379840,1056964608,3942645760,2801795072,67108864,1207959552,2617245696,3489660928,687865856,1694498816,2969567232,4244635648,2667577344,
3539992576,117440512,1258291200,2986344448,4261412864,704643072,1711276032,3305111552,2298478592,1560281088,285212672,3892314112,2751463424,1879048192,1006632960,721420288,
1728053248,3003121664,4278190080,100663296,1241513984,2650800128,3523215360,1895825408,1023410176,3909091328,2768240640,1543503872,268435456,3288334336,2281701376,3372220416,
2231369728,1358954496,486539264,3825205248,2818572288,2080374784,805306368,2466250752,3741319168,184549376,1191182336,3187671040,4060086272,637534208,1778384896,2097152000,
822083584,3841982464,2835349504,1342177280,469762048,3355443200,2214592512,654311424,1795162112,3204448256,4076863488,167772160,1174405120,2449473536,3724541952,3154116608,
4026531840,603979776,1744830464,2432696320,3707764736,150994944,1157627904,3858759680,2852126720,2113929216,838860800,3405774848,2264924160,1392508928,520093696,134217728,
1140850688,2415919104,3690987520,620756992,1761607680,3170893824,4043309056,1375731712,503316480,3388997632,2248146944,2130706432,855638016,3875536896,2868903936,587202560,
1862270976,3137339392,4143972352,234881024,1107296256,2516582400,3657433088,2030043136,889192448,3774873600,2902458368,1409286144,402653184,3422552064,2147483648,2533359616,
3674210304,251658240,1124073472,3120562176,4127195136,570425344,1845493760,3439329280,2164260864,1426063360,419430400,3758096384,2885681152,2013265920,872415232,1442840576,
436207616,3456106496,2181038080,2063597568,922746880,3808428032,2936012800,201326592,1073741824,2483027968,3623878656,553648128,1828716544,3103784960,4110417920,3791650816,
2919235584,2046820352,905969664,3472883712,2197815296,1459617792,452984832,3087007744,4093640704,536870912,1811939328,2499805184,3640655872,218103808,1090519040,0,
2399141888,50331648,2348810240,100663296,2298478592,83886080,2315255808,201326592,2197815296,251658240,2147483648,167772160,2231369728,150994944,2248146944,402653184,
2533359616,452984832,2483027968,503316480,2432696320,486539264,2449473536,335544320,2600468480,385875968,2550136832,301989888,2634022912,285212672,2650800128,805306368,
3204448256,855638016,3154116608,905969664,3103784960,889192448,3120562176,1006632960,3003121664,1056964608,2952790016,973078528,3036676096,956301312,3053453312,671088640,
2801795072,721420288,2751463424,771751936,2701131776,754974720,2717908992,603979776,2868903936,654311424,2818572288,570425344,2902458368,553648128,2919235584,1610612736,
4009754624,1660944384,3959422976,1711276032,3909091328,1694498816,3925868544,1811939328,3808428032,1862270976,3758096384,1778384896,3841982464,1761607680,3858759680,2013265920,
4143972352,2063597568,4093640704,2113929216,4043309056,2097152000,4060086272,1946157056,4211081216,1996488704,4160749568,1912602624,4244635648,1895825408,4261412864,1342177280,
3741319168,1392508928,3690987520,1442840576,3640655872,1426063360,3657433088,1543503872,3539992576,1593835520,3489660928,1509949440,3573547008,1493172224,3590324224,1207959552,
3338665984,1258291200,3288334336,1308622848,3238002688,1291845632,3254779904,1140850688,3405774848,1191182336,3355443200,1107296256,3439329280,1090519040,3456106496,3221225472,
1325400064,3271557120,1275068416,3321888768,1224736768,3305111552,1241513984,3422552064,1124073472,3472883712,1073741824,3388997632,1157627904,3372220416,1174405120,3623878656,
1459617792,3674210304,1409286144,3724541952,1358954496,3707764736,1375731712,3556769792,1526726656,3607101440,1476395008,3523215360,1560281088,3506438144,1577058304,4026531840,
2130706432,4076863488,2080374784,4127195136,2030043136,4110417920,2046820352,4227858432,1929379840,4278190080,1879048192,4194304000,1962934272,4177526784,1979711488,3892314112,
1728053248,3942645760,1677721600,3992977408,1627389952,3976200192,1644167168,3825205248,1795162112,3875536896,1744830464,3791650816,1828716544,3774873600,1845493760,2684354560,
788529152,2734686208,738197504,2785017856,687865856,2768240640,704643072,2885681152,587202560,2936012800,536870912,2852126720,620756992,2835349504,637534208,3087007744,
922746880,3137339392,872415232,3187671040,822083584,3170893824,838860800,3019898880,989855744,3070230528,939524096,2986344448,1023410176,2969567232,1040187392,2415919104,
520093696,2466250752,469762048,2516582400,419430400,2499805184,436207616,2617245696,318767104,2667577344,268435456,2583691264,352321536,2566914048,369098752,2281701376,
117440512,2332033024,67108864,2382364672,16777216,2365587456,33554432,2214592512,184549376,2264924160,134217728,2181038080,218103808,2164260864,234881024,0,
2634022912,654311424,3120562176,1308622848,3539992576,1761607680,4093640704,2617245696,16777216,3137339392,637534208,3523215360,1325400064,4110417920,1744830464,620756992,
3087007744,33554432,2667577344,1795162112,4127195136,1275068416,3506438144,3103784960,603979776,2650800128,50331648,4143972352,1778384896,3489660928,1291845632,1241513984,
3607101440,1828716544,4026531840,67108864,2566914048,587202560,3187671040,3590324224,1258291200,4043309056,1811939328,2550136832,83886080,3204448256,570425344,1862270976,
4060086272,1207959552,3573547008,553648128,3154116608,100663296,2600468480,4076863488,1845493760,3556769792,1224736768,3170893824,536870912,2583691264,117440512,2483027968,
150994944,3003121664,771751936,3657433088,1191182336,4244635648,1610612736,134217728,2499805184,788529152,2986344448,1174405120,3674210304,1627389952,4227858432,2969567232,
738197504,2516582400,184549376,4278190080,1644167168,3623878656,1157627904,754974720,2952790016,167772160,2533359616,1660944384,4261412864,1140850688,3640655872,3724541952,
1124073472,4177526784,1677721600,2415919104,218103808,3070230528,704643072,1107296256,3741319168,1694498816,4160749568,201326592,2432696320,721420288,3053453312,4211081216,
1711276032,3690987520,1090519040,3036676096,671088640,2449473536,251658240,1728053248,4194304000,1073741824,3707764736,687865856,3019898880,234881024,2466250752,889192448,
2818572288,301989888,2399141888,2063597568,3858759680,1543503872,3238002688,2835349504,872415232,2382364672,318767104,3875536896,2046820352,3221225472,1560281088,268435456,
2365587456,922746880,2852126720,1577058304,3271557120,2030043136,3825205248,2348810240,285212672,2868903936,905969664,3254779904,1593835520,3841982464,2013265920,2130706432,
3791650816,1476395008,3305111552,822083584,2885681152,369098752,2332033024,3808428032,2113929216,3288334336,1493172224,2902458368,805306368,2315255808,385875968,1509949440,
3338665984,2097152000,3758096384,335544320,2298478592,855638016,2919235584,3321888768,1526726656,3774873600,2080374784,2281701376,352321536,2936012800,838860800,2701131776,
1006632960,2248146944,452984832,4009754624,1912602624,3355443200,1426063360,1023410176,2684354560,436207616,2264924160,1929379840,3992977408,1409286144,3372220416,2214592512,
419430400,2734686208,1040187392,3388997632,1459617792,3976200192,1879048192,402653184,2231369728,1056964608,2717908992,1442840576,3405774848,1895825408,3959422976,3942645760,
1979711488,3422552064,1358954496,2768240640,939524096,2181038080,520093696,1996488704,3925868544,1342177280,3439329280,956301312,2751463424,503316480,2197815296,3456106496,
1392508928,3909091328,1946157056,2147483648,486539264,2801795072,973078528,1375731712,3472883712,1962934272,3892314112,469762048,2164260864,989855744,2785017856};

uint32_t ecc_buf[1];
uint32_t ecc_buf2[1];
unsigned int xi_tab[8] = {214,232,234,44,238,0,36,80};
unsigned int syn[2];
int cache[2];
struct bch_control ecc_bch={
.m = 8, 
.t = 1, 
.n = 255, 
.ecc_bytes = 1, 
.ecc_bits = 8, 
.a_pow_tab = a_pow_tab, 
.a_log_tab = a_log_tab, 
.mod8_tab = mod8_tab, 
.ecc_buf = ecc_buf, 
.ecc_buf2 = ecc_buf2, 
.xi_tab  = xi_tab,  
.syn = syn,  
.cache = cache, 
}; 
