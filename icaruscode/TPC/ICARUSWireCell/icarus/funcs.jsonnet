// This provides some util functions.

local g = import 'pgraph.jsonnet';

{
    //  Return a list of channel by anode index [1-8]
    // "AnodePlane:anode110"
    // "AnodePlane:anode120"
    // "AnodePlane:anode111"
    // "AnodePlane:anode121"
    // "AnodePlane:anode112"
    // "AnodePlane:anode122"
    // "AnodePlane:anode113"
    // "AnodePlane:anode123"
    local startch = [ [0, 2400, 8128],
                      [1152, 2400, 8128],
                      [13824, 16192, 21984],
                      [14976, 16192, 21984],
                      [27648, 30048, 35776],
                      [28800, 30048, 35776],
                      [41472, 43840, 49632],
                      [42624, 43840, 49632]
                    ],
    local wireplanes = [ 1055, 5599, 5599 ],
    anode_channels(n):: std.flattenArrays([std.range(startch[n][w], startch[n][w]+wireplanes[w]) for w in std.range(0,2)]),
    // anode_channels(n):: std.range(1056 * (n % 2) + 13312 * (n - n % 2) / 2, 1056 * (n % 2 + 1) - 1 + 13312 * (n - n % 2) / 2) + std.range(1056 * 2 + 13312 * (n - n % 2) / 2, 13312 - 1 + 13312 * (n - n % 2) / 2),

    // Return the number of split (1 or 2) for an anode
    anode_split(ident):: (ident%100 - ident%10)/10,

    //  Build a depofanout-[signal]-[framesummer]-[pipelines]-fanin graph.
    //  FrameSummer add up the two "split" anodes into one frame.
    //  Each branch of the pipelines operates on the summed signal frame.
    fansummer :: function(fout, sigpipes, summers, actpipes, fin, name="fansummer", outtags=[], tag_rules=[]) {

        local fanoutmult = std.length(sigpipes),
        local faninmult = std.length(actpipes),

        local fanout = g.pnode({
            type: fout,
            name: name,
            data: {
                multiplicity: fanoutmult,
                tag_rules: tag_rules,
            },
        }, nin=1, nout=fanoutmult),


        local fanin = g.pnode({
            type: fin,
            name: name,
            data: {
                multiplicity: faninmult,
                tags: outtags,
            },
        }, nin=faninmult, nout=1),

        local reducer = g.intern(innodes=sigpipes,
                                 outnodes=actpipes,
                                 centernodes=summers,
                                 edges= 
                                 // connecting signal and summer
                                 [g.edge(sigpipes[0], summers[0],0,0)]
                                 + [g.edge(sigpipes[1], summers[0],0,1)]
                                 + [g.edge(sigpipes[2], summers[1],0,0)]
                                 + [g.edge(sigpipes[3], summers[1],0,1)]
                                 + [g.edge(sigpipes[4], summers[2],0,0)]
                                 + [g.edge(sigpipes[5], summers[2],0,1)]
                                 + [g.edge(sigpipes[6], summers[3],0,0)]
                                 + [g.edge(sigpipes[7], summers[3],0,1)]
                                 // connecting summer and the operator pipelines
                                 + [g.edge(summers[n], actpipes[n]) for n in std.range(0,faninmult-1)],
                                 name=name),

        ret: g.intern(innodes=[fanout],
                      outnodes=[fanin],
                      centernodes=[reducer],
                      edges=
                      [g.edge(fanout, sigpipes[n], n, 0) for n in std.range(0, fanoutmult-1)] +
                      [g.edge(reducer, fanin, n, n) for n in std.range(0, faninmult-1)],
                      name=name),
    }.ret,

    //  Build a depofanout-[signal]-[framesummer]-[pipelines]-fanin graph.
    //  FrameSummer add up the two "split" anodes into one frame.
    //  Each branch of the pipelines operates on the summed signal frame.
    fansummeryz :: function(fout, sigpipes, summers, actpipes, fin, name="fansummeryz", outtags=[], tag_rules=[]) {

        local fanoutmult = std.length(sigpipes),
        local faninmult = std.length(actpipes),

        local fanout = g.pnode({
            type: fout,
            name: name,
            data: {
                multiplicity: fanoutmult,
                tag_rules: tag_rules,
            },
        }, nin=1, nout=fanoutmult),


        local fanin = g.pnode({
            type: fin,
            name: name,
            data: {
                multiplicity: faninmult,
                tags: outtags,
            },
        }, nin=faninmult, nout=1),

        local reduceryz = g.intern(innodes=sigpipes,
                                 outnodes=actpipes,
                                 centernodes=summers,
                                 edges= 
                                 // connecting signal and summer
                                 [g.edge(sigpipes[0], summers[0],0,0)]
                                   + [g.edge(sigpipes[1], summers[0],0,1)]
                                   + [g.edge(sigpipes[2], summers[0],0,2)]
                                   + [g.edge(sigpipes[3], summers[0],0,3)]
                                   + [g.edge(sigpipes[4], summers[0],0,4)]
                                   + [g.edge(sigpipes[5], summers[0],0,5)]
                                   + [g.edge(sigpipes[6], summers[0],0,6)]
                                   + [g.edge(sigpipes[7], summers[0],0,7)]
                                   + [g.edge(sigpipes[8], summers[0],0,8)]
                                   + [g.edge(sigpipes[9], summers[0],0,9)]
                                   + [g.edge(sigpipes[10], summers[0],0,10)]
                                   + [g.edge(sigpipes[11], summers[0],0,11)]
                                   + [g.edge(sigpipes[12], summers[0],0,12)]
                                   + [g.edge(sigpipes[13], summers[0],0,13)]
                                   + [g.edge(sigpipes[14], summers[0],0,14)]
                                   + [g.edge(sigpipes[15], summers[0],0,15)]
                                   + [g.edge(sigpipes[16], summers[0],0,16)]
                                   + [g.edge(sigpipes[17], summers[0],0,17)]
                                   + [g.edge(sigpipes[18], summers[0],0,18)]
                                   + [g.edge(sigpipes[19], summers[0],0,19)]
                                   + [g.edge(sigpipes[20], summers[0],0,20)]
                                   + [g.edge(sigpipes[21], summers[0],0,21)]
                                   + [g.edge(sigpipes[22], summers[0],0,22)]
                                   + [g.edge(sigpipes[23], summers[0],0,23)]
                                   + [g.edge(sigpipes[24], summers[0],0,24)]
                                   + [g.edge(sigpipes[25], summers[0],0,25)]
                                   + [g.edge(sigpipes[26], summers[0],0,26)]
                                   + [g.edge(sigpipes[27], summers[0],0,27)]
                                   + [g.edge(sigpipes[28], summers[0],0,28)]
                                   + [g.edge(sigpipes[29], summers[0],0,29)]
                                   + [g.edge(sigpipes[30], summers[0],0,30)]
                                   + [g.edge(sigpipes[31], summers[0],0,31)]
                                   + [g.edge(sigpipes[32], summers[0],0,32)]
                                   + [g.edge(sigpipes[33], summers[0],0,33)]
                                   + [g.edge(sigpipes[34], summers[0],0,34)]
                                   + [g.edge(sigpipes[35], summers[0],0,35)]
                                   + [g.edge(sigpipes[36], summers[0],0,36)]
                                   + [g.edge(sigpipes[37], summers[0],0,37)]
                                   + [g.edge(sigpipes[38], summers[0],0,38)]
                                   + [g.edge(sigpipes[39], summers[0],0,39)]
                                   + [g.edge(sigpipes[40], summers[0],0,40)]
                                   + [g.edge(sigpipes[41], summers[0],0,41)]
                                   + [g.edge(sigpipes[42], summers[0],0,42)]
                                   + [g.edge(sigpipes[43], summers[0],0,43)]
                                   + [g.edge(sigpipes[44], summers[0],0,44)]
                                   + [g.edge(sigpipes[45], summers[0],0,45)]
                                   + [g.edge(sigpipes[46], summers[0],0,46)]
                                   + [g.edge(sigpipes[47], summers[0],0,47)]
                                   + [g.edge(sigpipes[48], summers[0],0,48)]
                                   + [g.edge(sigpipes[49], summers[0],0,49)]
                                   + [g.edge(sigpipes[50], summers[0],0,50)]
                                   + [g.edge(sigpipes[51], summers[0],0,51)]
                                   + [g.edge(sigpipes[52], summers[0],0,52)]
                                   + [g.edge(sigpipes[53], summers[0],0,53)]
                                   + [g.edge(sigpipes[54], summers[0],0,54)]
                                   + [g.edge(sigpipes[55], summers[0],0,55)]
                                   + [g.edge(sigpipes[56], summers[0],0,56)]
                                   + [g.edge(sigpipes[57], summers[0],0,57)]
                                   + [g.edge(sigpipes[58], summers[0],0,58)]
                                   + [g.edge(sigpipes[59], summers[0],0,59)]
                                   + [g.edge(sigpipes[60], summers[0],0,60)]
                                   + [g.edge(sigpipes[61], summers[0],0,61)]
                                   + [g.edge(sigpipes[62], summers[0],0,62)]
                                   + [g.edge(sigpipes[63], summers[0],0,63)]
                                   + [g.edge(sigpipes[64], summers[0],0,64)]
                                   + [g.edge(sigpipes[65], summers[0],0,65)]
                                   + [g.edge(sigpipes[66], summers[0],0,66)]
                                   + [g.edge(sigpipes[67], summers[0],0,67)]
                                   + [g.edge(sigpipes[68], summers[0],0,68)]
                                   + [g.edge(sigpipes[69], summers[0],0,69)]
                                   + [g.edge(sigpipes[70], summers[0],0,70)]
                                   + [g.edge(sigpipes[71], summers[0],0,71)]
                                   + [g.edge(sigpipes[72], summers[0],0,72)]
                                   + [g.edge(sigpipes[73], summers[0],0,73)]
                                   + [g.edge(sigpipes[74], summers[0],0,74)]
                                   + [g.edge(sigpipes[75], summers[0],0,75)]
                                   + [g.edge(sigpipes[76], summers[0],0,76)]
                                   + [g.edge(sigpipes[77], summers[0],0,77)]
                                   + [g.edge(sigpipes[78], summers[0],0,78)]
                                   + [g.edge(sigpipes[79], summers[0],0,79)]
                                   + [g.edge(sigpipes[80], summers[0],0,80)]
                                   + [g.edge(sigpipes[81], summers[0],0,81)]
                                   + [g.edge(sigpipes[82], summers[0],0,82)]
                                   + [g.edge(sigpipes[83], summers[0],0,83)]
                                   + [g.edge(sigpipes[84], summers[0],0,84)]
                                   + [g.edge(sigpipes[85], summers[0],0,85)]
                                   + [g.edge(sigpipes[86], summers[0],0,86)]
                                   + [g.edge(sigpipes[87], summers[0],0,87)]
                                   + [g.edge(sigpipes[88], summers[0],0,88)]
                                   + [g.edge(sigpipes[89], summers[0],0,89)]
                                   + [g.edge(sigpipes[90], summers[1],0,0)]
                                   + [g.edge(sigpipes[91], summers[1],0,1)]
                                   + [g.edge(sigpipes[92], summers[1],0,2)]
                                   + [g.edge(sigpipes[93], summers[1],0,3)]
                                   + [g.edge(sigpipes[94], summers[1],0,4)]
                                   + [g.edge(sigpipes[95], summers[1],0,5)]
                                   + [g.edge(sigpipes[96], summers[1],0,6)]
                                   + [g.edge(sigpipes[97], summers[1],0,7)]
                                   + [g.edge(sigpipes[98], summers[1],0,8)]
                                   + [g.edge(sigpipes[99], summers[1],0,9)]
                                   + [g.edge(sigpipes[100], summers[1],0,10)]
                                   + [g.edge(sigpipes[101], summers[1],0,11)]
                                   + [g.edge(sigpipes[102], summers[1],0,12)]
                                   + [g.edge(sigpipes[103], summers[1],0,13)]
                                   + [g.edge(sigpipes[104], summers[1],0,14)]
                                   + [g.edge(sigpipes[105], summers[1],0,15)]
                                   + [g.edge(sigpipes[106], summers[1],0,16)]
                                   + [g.edge(sigpipes[107], summers[1],0,17)]
                                   + [g.edge(sigpipes[108], summers[1],0,18)]
                                   + [g.edge(sigpipes[109], summers[1],0,19)]
                                   + [g.edge(sigpipes[110], summers[1],0,20)]
                                   + [g.edge(sigpipes[111], summers[1],0,21)]
                                   + [g.edge(sigpipes[112], summers[1],0,22)]
                                   + [g.edge(sigpipes[113], summers[1],0,23)]
                                   + [g.edge(sigpipes[114], summers[1],0,24)]
                                   + [g.edge(sigpipes[115], summers[1],0,25)]
                                   + [g.edge(sigpipes[116], summers[1],0,26)]
                                   + [g.edge(sigpipes[117], summers[1],0,27)]
                                   + [g.edge(sigpipes[118], summers[1],0,28)]
                                   + [g.edge(sigpipes[119], summers[1],0,29)]
                                   + [g.edge(sigpipes[120], summers[1],0,30)]
                                   + [g.edge(sigpipes[121], summers[1],0,31)]
                                   + [g.edge(sigpipes[122], summers[1],0,32)]
                                   + [g.edge(sigpipes[123], summers[1],0,33)]
                                   + [g.edge(sigpipes[124], summers[1],0,34)]
                                   + [g.edge(sigpipes[125], summers[1],0,35)]
                                   + [g.edge(sigpipes[126], summers[1],0,36)]
                                   + [g.edge(sigpipes[127], summers[1],0,37)]
                                   + [g.edge(sigpipes[128], summers[1],0,38)]
                                   + [g.edge(sigpipes[129], summers[1],0,39)]
                                   + [g.edge(sigpipes[130], summers[1],0,40)]
                                   + [g.edge(sigpipes[131], summers[1],0,41)]
                                   + [g.edge(sigpipes[132], summers[1],0,42)]
                                   + [g.edge(sigpipes[133], summers[1],0,43)]
                                   + [g.edge(sigpipes[134], summers[1],0,44)]
                                   + [g.edge(sigpipes[135], summers[1],0,45)]
                                   + [g.edge(sigpipes[136], summers[1],0,46)]
                                   + [g.edge(sigpipes[137], summers[1],0,47)]
                                   + [g.edge(sigpipes[138], summers[1],0,48)]
                                   + [g.edge(sigpipes[139], summers[1],0,49)]
                                   + [g.edge(sigpipes[140], summers[1],0,50)]
                                   + [g.edge(sigpipes[141], summers[1],0,51)]
                                   + [g.edge(sigpipes[142], summers[1],0,52)]
                                   + [g.edge(sigpipes[143], summers[1],0,53)]
                                   + [g.edge(sigpipes[144], summers[1],0,54)]
                                   + [g.edge(sigpipes[145], summers[1],0,55)]
                                   + [g.edge(sigpipes[146], summers[1],0,56)]
                                   + [g.edge(sigpipes[147], summers[1],0,57)]
                                   + [g.edge(sigpipes[148], summers[1],0,58)]
                                   + [g.edge(sigpipes[149], summers[1],0,59)]
                                   + [g.edge(sigpipes[150], summers[1],0,60)]
                                   + [g.edge(sigpipes[151], summers[1],0,61)]
                                   + [g.edge(sigpipes[152], summers[1],0,62)]
                                   + [g.edge(sigpipes[153], summers[1],0,63)]
                                   + [g.edge(sigpipes[154], summers[1],0,64)]
                                   + [g.edge(sigpipes[155], summers[1],0,65)]
                                   + [g.edge(sigpipes[156], summers[1],0,66)]
                                   + [g.edge(sigpipes[157], summers[1],0,67)]
                                   + [g.edge(sigpipes[158], summers[1],0,68)]
                                   + [g.edge(sigpipes[159], summers[1],0,69)]
                                   + [g.edge(sigpipes[160], summers[1],0,70)]
                                   + [g.edge(sigpipes[161], summers[1],0,71)]
                                   + [g.edge(sigpipes[162], summers[1],0,72)]
                                   + [g.edge(sigpipes[163], summers[1],0,73)]
                                   + [g.edge(sigpipes[164], summers[1],0,74)]
                                   + [g.edge(sigpipes[165], summers[1],0,75)]
                                   + [g.edge(sigpipes[166], summers[1],0,76)]
                                   + [g.edge(sigpipes[167], summers[1],0,77)]
                                   + [g.edge(sigpipes[168], summers[1],0,78)]
                                   + [g.edge(sigpipes[169], summers[1],0,79)]
                                   + [g.edge(sigpipes[170], summers[1],0,80)]
                                   + [g.edge(sigpipes[171], summers[1],0,81)]
                                   + [g.edge(sigpipes[172], summers[1],0,82)]
                                   + [g.edge(sigpipes[173], summers[1],0,83)]
                                   + [g.edge(sigpipes[174], summers[1],0,84)]
                                   + [g.edge(sigpipes[175], summers[1],0,85)]
                                   + [g.edge(sigpipes[176], summers[1],0,86)]
                                   + [g.edge(sigpipes[177], summers[1],0,87)]
                                   + [g.edge(sigpipes[178], summers[1],0,88)]
                                   + [g.edge(sigpipes[179], summers[1],0,89)]
                                   + [g.edge(sigpipes[180], summers[2],0,0)]
                                   + [g.edge(sigpipes[181], summers[2],0,1)]
                                   + [g.edge(sigpipes[182], summers[2],0,2)]
                                   + [g.edge(sigpipes[183], summers[2],0,3)]
                                   + [g.edge(sigpipes[184], summers[2],0,4)]
                                   + [g.edge(sigpipes[185], summers[2],0,5)]
                                   + [g.edge(sigpipes[186], summers[2],0,6)]
                                   + [g.edge(sigpipes[187], summers[2],0,7)]
                                   + [g.edge(sigpipes[188], summers[2],0,8)]
                                   + [g.edge(sigpipes[189], summers[2],0,9)]
                                   + [g.edge(sigpipes[190], summers[2],0,10)]
                                   + [g.edge(sigpipes[191], summers[2],0,11)]
                                   + [g.edge(sigpipes[192], summers[2],0,12)]
                                   + [g.edge(sigpipes[193], summers[2],0,13)]
                                   + [g.edge(sigpipes[194], summers[2],0,14)]
                                   + [g.edge(sigpipes[195], summers[2],0,15)]
                                   + [g.edge(sigpipes[196], summers[2],0,16)]
                                   + [g.edge(sigpipes[197], summers[2],0,17)]
                                   + [g.edge(sigpipes[198], summers[2],0,18)]
                                   + [g.edge(sigpipes[199], summers[2],0,19)]
                                   + [g.edge(sigpipes[200], summers[2],0,20)]
                                   + [g.edge(sigpipes[201], summers[2],0,21)]
                                   + [g.edge(sigpipes[202], summers[2],0,22)]
                                   + [g.edge(sigpipes[203], summers[2],0,23)]
                                   + [g.edge(sigpipes[204], summers[2],0,24)]
                                   + [g.edge(sigpipes[205], summers[2],0,25)]
                                   + [g.edge(sigpipes[206], summers[2],0,26)]
                                   + [g.edge(sigpipes[207], summers[2],0,27)]
                                   + [g.edge(sigpipes[208], summers[2],0,28)]
                                   + [g.edge(sigpipes[209], summers[2],0,29)]
                                   + [g.edge(sigpipes[210], summers[2],0,30)]
                                   + [g.edge(sigpipes[211], summers[2],0,31)]
                                   + [g.edge(sigpipes[212], summers[2],0,32)]
                                   + [g.edge(sigpipes[213], summers[2],0,33)]
                                   + [g.edge(sigpipes[214], summers[2],0,34)]
                                   + [g.edge(sigpipes[215], summers[2],0,35)]
                                   + [g.edge(sigpipes[216], summers[2],0,36)]
                                   + [g.edge(sigpipes[217], summers[2],0,37)]
                                   + [g.edge(sigpipes[218], summers[2],0,38)]
                                   + [g.edge(sigpipes[219], summers[2],0,39)]
                                   + [g.edge(sigpipes[220], summers[2],0,40)]
                                   + [g.edge(sigpipes[221], summers[2],0,41)]
                                   + [g.edge(sigpipes[222], summers[2],0,42)]
                                   + [g.edge(sigpipes[223], summers[2],0,43)]
                                   + [g.edge(sigpipes[224], summers[2],0,44)]
                                   + [g.edge(sigpipes[225], summers[2],0,45)]
                                   + [g.edge(sigpipes[226], summers[2],0,46)]
                                   + [g.edge(sigpipes[227], summers[2],0,47)]
                                   + [g.edge(sigpipes[228], summers[2],0,48)]
                                   + [g.edge(sigpipes[229], summers[2],0,49)]
                                   + [g.edge(sigpipes[230], summers[2],0,50)]
                                   + [g.edge(sigpipes[231], summers[2],0,51)]
                                   + [g.edge(sigpipes[232], summers[2],0,52)]
                                   + [g.edge(sigpipes[233], summers[2],0,53)]
                                   + [g.edge(sigpipes[234], summers[2],0,54)]
                                   + [g.edge(sigpipes[235], summers[2],0,55)]
                                   + [g.edge(sigpipes[236], summers[2],0,56)]
                                   + [g.edge(sigpipes[237], summers[2],0,57)]
                                   + [g.edge(sigpipes[238], summers[2],0,58)]
                                   + [g.edge(sigpipes[239], summers[2],0,59)]
                                   + [g.edge(sigpipes[240], summers[2],0,60)]
                                   + [g.edge(sigpipes[241], summers[2],0,61)]
                                   + [g.edge(sigpipes[242], summers[2],0,62)]
                                   + [g.edge(sigpipes[243], summers[2],0,63)]
                                   + [g.edge(sigpipes[244], summers[2],0,64)]
                                   + [g.edge(sigpipes[245], summers[2],0,65)]
                                   + [g.edge(sigpipes[246], summers[2],0,66)]
                                   + [g.edge(sigpipes[247], summers[2],0,67)]
                                   + [g.edge(sigpipes[248], summers[2],0,68)]
                                   + [g.edge(sigpipes[249], summers[2],0,69)]
                                   + [g.edge(sigpipes[250], summers[2],0,70)]
                                   + [g.edge(sigpipes[251], summers[2],0,71)]
                                   + [g.edge(sigpipes[252], summers[2],0,72)]
                                   + [g.edge(sigpipes[253], summers[2],0,73)]
                                   + [g.edge(sigpipes[254], summers[2],0,74)]
                                   + [g.edge(sigpipes[255], summers[2],0,75)]
                                   + [g.edge(sigpipes[256], summers[2],0,76)]
                                   + [g.edge(sigpipes[257], summers[2],0,77)]
                                   + [g.edge(sigpipes[258], summers[2],0,78)]
                                   + [g.edge(sigpipes[259], summers[2],0,79)]
                                   + [g.edge(sigpipes[260], summers[2],0,80)]
                                   + [g.edge(sigpipes[261], summers[2],0,81)]
                                   + [g.edge(sigpipes[262], summers[2],0,82)]
                                   + [g.edge(sigpipes[263], summers[2],0,83)]
                                   + [g.edge(sigpipes[264], summers[2],0,84)]
                                   + [g.edge(sigpipes[265], summers[2],0,85)]
                                   + [g.edge(sigpipes[266], summers[2],0,86)]
                                   + [g.edge(sigpipes[267], summers[2],0,87)]
                                   + [g.edge(sigpipes[268], summers[2],0,88)]
                                   + [g.edge(sigpipes[269], summers[2],0,89)]
                                   + [g.edge(sigpipes[270], summers[3],0,0)]
                                   + [g.edge(sigpipes[271], summers[3],0,1)]
                                   + [g.edge(sigpipes[272], summers[3],0,2)]
                                   + [g.edge(sigpipes[273], summers[3],0,3)]
                                   + [g.edge(sigpipes[274], summers[3],0,4)]
                                   + [g.edge(sigpipes[275], summers[3],0,5)]
                                   + [g.edge(sigpipes[276], summers[3],0,6)]
                                   + [g.edge(sigpipes[277], summers[3],0,7)]
                                   + [g.edge(sigpipes[278], summers[3],0,8)]
                                   + [g.edge(sigpipes[279], summers[3],0,9)]
                                   + [g.edge(sigpipes[280], summers[3],0,10)]
                                   + [g.edge(sigpipes[281], summers[3],0,11)]
                                   + [g.edge(sigpipes[282], summers[3],0,12)]
                                   + [g.edge(sigpipes[283], summers[3],0,13)]
                                   + [g.edge(sigpipes[284], summers[3],0,14)]
                                   + [g.edge(sigpipes[285], summers[3],0,15)]
                                   + [g.edge(sigpipes[286], summers[3],0,16)]
                                   + [g.edge(sigpipes[287], summers[3],0,17)]
                                   + [g.edge(sigpipes[288], summers[3],0,18)]
                                   + [g.edge(sigpipes[289], summers[3],0,19)]
                                   + [g.edge(sigpipes[290], summers[3],0,20)]
                                   + [g.edge(sigpipes[291], summers[3],0,21)]
                                   + [g.edge(sigpipes[292], summers[3],0,22)]
                                   + [g.edge(sigpipes[293], summers[3],0,23)]
                                   + [g.edge(sigpipes[294], summers[3],0,24)]
                                   + [g.edge(sigpipes[295], summers[3],0,25)]
                                   + [g.edge(sigpipes[296], summers[3],0,26)]
                                   + [g.edge(sigpipes[297], summers[3],0,27)]
                                   + [g.edge(sigpipes[298], summers[3],0,28)]
                                   + [g.edge(sigpipes[299], summers[3],0,29)]
                                   + [g.edge(sigpipes[300], summers[3],0,30)]
                                   + [g.edge(sigpipes[301], summers[3],0,31)]
                                   + [g.edge(sigpipes[302], summers[3],0,32)]
                                   + [g.edge(sigpipes[303], summers[3],0,33)]
                                   + [g.edge(sigpipes[304], summers[3],0,34)]
                                   + [g.edge(sigpipes[305], summers[3],0,35)]
                                   + [g.edge(sigpipes[306], summers[3],0,36)]
                                   + [g.edge(sigpipes[307], summers[3],0,37)]
                                   + [g.edge(sigpipes[308], summers[3],0,38)]
                                   + [g.edge(sigpipes[309], summers[3],0,39)]
                                   + [g.edge(sigpipes[310], summers[3],0,40)]
                                   + [g.edge(sigpipes[311], summers[3],0,41)]
                                   + [g.edge(sigpipes[312], summers[3],0,42)]
                                   + [g.edge(sigpipes[313], summers[3],0,43)]
                                   + [g.edge(sigpipes[314], summers[3],0,44)]
                                   + [g.edge(sigpipes[315], summers[3],0,45)]
                                   + [g.edge(sigpipes[316], summers[3],0,46)]
                                   + [g.edge(sigpipes[317], summers[3],0,47)]
                                   + [g.edge(sigpipes[318], summers[3],0,48)]
                                   + [g.edge(sigpipes[319], summers[3],0,49)]
                                   + [g.edge(sigpipes[320], summers[3],0,50)]
                                   + [g.edge(sigpipes[321], summers[3],0,51)]
                                   + [g.edge(sigpipes[322], summers[3],0,52)]
                                   + [g.edge(sigpipes[323], summers[3],0,53)]
                                   + [g.edge(sigpipes[324], summers[3],0,54)]
                                   + [g.edge(sigpipes[325], summers[3],0,55)]
                                   + [g.edge(sigpipes[326], summers[3],0,56)]
                                   + [g.edge(sigpipes[327], summers[3],0,57)]
                                   + [g.edge(sigpipes[328], summers[3],0,58)]
                                   + [g.edge(sigpipes[329], summers[3],0,59)]
                                   + [g.edge(sigpipes[330], summers[3],0,60)]
                                   + [g.edge(sigpipes[331], summers[3],0,61)]
                                   + [g.edge(sigpipes[332], summers[3],0,62)]
                                   + [g.edge(sigpipes[333], summers[3],0,63)]
                                   + [g.edge(sigpipes[334], summers[3],0,64)]
                                   + [g.edge(sigpipes[335], summers[3],0,65)]
                                   + [g.edge(sigpipes[336], summers[3],0,66)]
                                   + [g.edge(sigpipes[337], summers[3],0,67)]
                                   + [g.edge(sigpipes[338], summers[3],0,68)]
                                   + [g.edge(sigpipes[339], summers[3],0,69)]
                                   + [g.edge(sigpipes[340], summers[3],0,70)]
                                   + [g.edge(sigpipes[341], summers[3],0,71)]
                                   + [g.edge(sigpipes[342], summers[3],0,72)]
                                   + [g.edge(sigpipes[343], summers[3],0,73)]
                                   + [g.edge(sigpipes[344], summers[3],0,74)]
                                   + [g.edge(sigpipes[345], summers[3],0,75)]
                                   + [g.edge(sigpipes[346], summers[3],0,76)]
                                   + [g.edge(sigpipes[347], summers[3],0,77)]
                                   + [g.edge(sigpipes[348], summers[3],0,78)]
                                   + [g.edge(sigpipes[349], summers[3],0,79)]
                                   + [g.edge(sigpipes[350], summers[3],0,80)]
                                   + [g.edge(sigpipes[351], summers[3],0,81)]
                                   + [g.edge(sigpipes[352], summers[3],0,82)]
                                   + [g.edge(sigpipes[353], summers[3],0,83)]
                                   + [g.edge(sigpipes[354], summers[3],0,84)]
                                   + [g.edge(sigpipes[355], summers[3],0,85)]
                                   + [g.edge(sigpipes[356], summers[3],0,86)]
                                   + [g.edge(sigpipes[357], summers[3],0,87)]
                                   + [g.edge(sigpipes[358], summers[3],0,88)]
                                   + [g.edge(sigpipes[359], summers[3],0,89)]
                                 // connecting summer and the operator pipelines
                                 + [g.edge(summers[n], actpipes[n]) for n in std.range(0,faninmult-1)],
                                 name=name),

        ret: g.intern(innodes=[fanout],
                      outnodes=[fanin],
                      centernodes=[reduceryz],
                      edges=
                      [g.edge(fanout, sigpipes[n], n, 0) for n in std.range(0, fanoutmult-1)] +
                      [g.edge(reduceryz, fanin, n, n) for n in std.range(0, faninmult-1)],
                      name=name),
    }.ret,



     //  Build a depofanout-[drift]-[signal]-[framesummer]-[pipelines]-fanin graph.
    //  FrameSummer add up the two "split" anodes into one frame.
    //  Each branch of the pipelines operates on the summed signal frame.
    fandrifter :: function(fout,driftpipes, sigpipes, summers, actpipes, fin, name="fandrifter", outtags=[], tag_rules=[]) {

        local fanoutmult = std.length(driftpipes),
        local faninmult = std.length(actpipes),

        local fanout = g.pnode({
            type: fout,
            name: name,
            data: {
                multiplicity: fanoutmult,
                tag_rules: tag_rules,
            },
        }, nin=1, nout=fanoutmult),


        local fanin = g.pnode({
            type: fin,
            name: name,
            data: {
                multiplicity: faninmult,
                tags: outtags,
            },
        }, nin=faninmult, nout=1),

	local drift = g.intern(innodes=driftpipes,
                                 outnodes=actpipes,
                                 centernodes=sigpipes+summers,
                                 edges=
				 [g.edge(driftpipes[n], sigpipes[n]) for n in std.range(0,fanoutmult-1)]
                                 // connecting signal and summer
                                   + [g.edge(sigpipes[0], summers[0],0,0)]
                                   + [g.edge(sigpipes[1], summers[0],0,1)]
                                   + [g.edge(sigpipes[2], summers[0],0,2)]
                                   + [g.edge(sigpipes[3], summers[0],0,3)]
                                   + [g.edge(sigpipes[4], summers[0],0,4)]
                                   + [g.edge(sigpipes[5], summers[0],0,5)]
                                   + [g.edge(sigpipes[6], summers[0],0,6)]
                                   + [g.edge(sigpipes[7], summers[0],0,7)]
                                   + [g.edge(sigpipes[8], summers[0],0,8)]
                                   + [g.edge(sigpipes[9], summers[0],0,9)]
                                   + [g.edge(sigpipes[10], summers[0],0,10)]
                                   + [g.edge(sigpipes[11], summers[0],0,11)]
                                   + [g.edge(sigpipes[12], summers[0],0,12)]
                                   + [g.edge(sigpipes[13], summers[0],0,13)]
                                   + [g.edge(sigpipes[14], summers[0],0,14)]
                                   + [g.edge(sigpipes[15], summers[0],0,15)]
                                   + [g.edge(sigpipes[16], summers[0],0,16)]
                                   + [g.edge(sigpipes[17], summers[0],0,17)]
                                   + [g.edge(sigpipes[18], summers[0],0,18)]
                                   + [g.edge(sigpipes[19], summers[0],0,19)]
                                   + [g.edge(sigpipes[20], summers[0],0,20)]
                                   + [g.edge(sigpipes[21], summers[0],0,21)]
                                   + [g.edge(sigpipes[22], summers[0],0,22)]
                                   + [g.edge(sigpipes[23], summers[0],0,23)]
                                   + [g.edge(sigpipes[24], summers[0],0,24)]
                                   + [g.edge(sigpipes[25], summers[0],0,25)]
                                   + [g.edge(sigpipes[26], summers[0],0,26)]
                                   + [g.edge(sigpipes[27], summers[0],0,27)]
                                   + [g.edge(sigpipes[28], summers[0],0,28)]
                                   + [g.edge(sigpipes[29], summers[0],0,29)]
                                   + [g.edge(sigpipes[30], summers[0],0,30)]
                                   + [g.edge(sigpipes[31], summers[0],0,31)]
                                   + [g.edge(sigpipes[32], summers[0],0,32)]
                                   + [g.edge(sigpipes[33], summers[0],0,33)]
                                   + [g.edge(sigpipes[34], summers[0],0,34)]
                                   + [g.edge(sigpipes[35], summers[0],0,35)]
                                   + [g.edge(sigpipes[36], summers[0],0,36)]
                                   + [g.edge(sigpipes[37], summers[0],0,37)]
                                   + [g.edge(sigpipes[38], summers[0],0,38)]
                                   + [g.edge(sigpipes[39], summers[0],0,39)]
                                   + [g.edge(sigpipes[40], summers[0],0,40)]
                                   + [g.edge(sigpipes[41], summers[0],0,41)]
                                   + [g.edge(sigpipes[42], summers[0],0,42)]
                                   + [g.edge(sigpipes[43], summers[0],0,43)]
                                   + [g.edge(sigpipes[44], summers[0],0,44)]
                                   + [g.edge(sigpipes[45], summers[0],0,45)]
                                   + [g.edge(sigpipes[46], summers[0],0,46)]
                                   + [g.edge(sigpipes[47], summers[0],0,47)]
                                   + [g.edge(sigpipes[48], summers[0],0,48)]
                                   + [g.edge(sigpipes[49], summers[0],0,49)]
                                   + [g.edge(sigpipes[50], summers[0],0,50)]
                                   + [g.edge(sigpipes[51], summers[0],0,51)]
                                   + [g.edge(sigpipes[52], summers[0],0,52)]
                                   + [g.edge(sigpipes[53], summers[0],0,53)]
                                   + [g.edge(sigpipes[54], summers[0],0,54)]
                                   + [g.edge(sigpipes[55], summers[0],0,55)]
                                   + [g.edge(sigpipes[56], summers[0],0,56)]
                                   + [g.edge(sigpipes[57], summers[0],0,57)]
                                   + [g.edge(sigpipes[58], summers[0],0,58)]
                                   + [g.edge(sigpipes[59], summers[0],0,59)]
                                   + [g.edge(sigpipes[60], summers[0],0,60)]
                                   + [g.edge(sigpipes[61], summers[0],0,61)]
                                   + [g.edge(sigpipes[62], summers[0],0,62)]
                                   + [g.edge(sigpipes[63], summers[0],0,63)]
                                   + [g.edge(sigpipes[64], summers[0],0,64)]
                                   + [g.edge(sigpipes[65], summers[0],0,65)]
                                   + [g.edge(sigpipes[66], summers[0],0,66)]
                                   + [g.edge(sigpipes[67], summers[0],0,67)]
                                   + [g.edge(sigpipes[68], summers[0],0,68)]
                                   + [g.edge(sigpipes[69], summers[0],0,69)]
                                   + [g.edge(sigpipes[70], summers[0],0,70)]
                                   + [g.edge(sigpipes[71], summers[0],0,71)]
                                   + [g.edge(sigpipes[72], summers[0],0,72)]
                                   + [g.edge(sigpipes[73], summers[0],0,73)]
                                   + [g.edge(sigpipes[74], summers[0],0,74)]
                                   + [g.edge(sigpipes[75], summers[0],0,75)]
                                   + [g.edge(sigpipes[76], summers[0],0,76)]
                                   + [g.edge(sigpipes[77], summers[0],0,77)]
                                   + [g.edge(sigpipes[78], summers[0],0,78)]
                                   + [g.edge(sigpipes[79], summers[0],0,79)]
                                   + [g.edge(sigpipes[80], summers[0],0,80)]
                                   + [g.edge(sigpipes[81], summers[0],0,81)]
                                   + [g.edge(sigpipes[82], summers[0],0,82)]
                                   + [g.edge(sigpipes[83], summers[0],0,83)]
                                   + [g.edge(sigpipes[84], summers[0],0,84)]
                                   + [g.edge(sigpipes[85], summers[0],0,85)]
                                   + [g.edge(sigpipes[86], summers[0],0,86)]
                                   + [g.edge(sigpipes[87], summers[0],0,87)]
                                   + [g.edge(sigpipes[88], summers[0],0,88)]
                                   + [g.edge(sigpipes[89], summers[0],0,89)]
                                   + [g.edge(sigpipes[90], summers[1],0,0)]
                                   + [g.edge(sigpipes[91], summers[1],0,1)]
                                   + [g.edge(sigpipes[92], summers[1],0,2)]
                                   + [g.edge(sigpipes[93], summers[1],0,3)]
                                   + [g.edge(sigpipes[94], summers[1],0,4)]
                                   + [g.edge(sigpipes[95], summers[1],0,5)]
                                   + [g.edge(sigpipes[96], summers[1],0,6)]
                                   + [g.edge(sigpipes[97], summers[1],0,7)]
                                   + [g.edge(sigpipes[98], summers[1],0,8)]
                                   + [g.edge(sigpipes[99], summers[1],0,9)]
                                   + [g.edge(sigpipes[100], summers[1],0,10)]
                                   + [g.edge(sigpipes[101], summers[1],0,11)]
                                   + [g.edge(sigpipes[102], summers[1],0,12)]
                                   + [g.edge(sigpipes[103], summers[1],0,13)]
                                   + [g.edge(sigpipes[104], summers[1],0,14)]
                                   + [g.edge(sigpipes[105], summers[1],0,15)]
                                   + [g.edge(sigpipes[106], summers[1],0,16)]
                                   + [g.edge(sigpipes[107], summers[1],0,17)]
                                   + [g.edge(sigpipes[108], summers[1],0,18)]
                                   + [g.edge(sigpipes[109], summers[1],0,19)]
                                   + [g.edge(sigpipes[110], summers[1],0,20)]
                                   + [g.edge(sigpipes[111], summers[1],0,21)]
                                   + [g.edge(sigpipes[112], summers[1],0,22)]
                                   + [g.edge(sigpipes[113], summers[1],0,23)]
                                   + [g.edge(sigpipes[114], summers[1],0,24)]
                                   + [g.edge(sigpipes[115], summers[1],0,25)]
                                   + [g.edge(sigpipes[116], summers[1],0,26)]
                                   + [g.edge(sigpipes[117], summers[1],0,27)]
                                   + [g.edge(sigpipes[118], summers[1],0,28)]
                                   + [g.edge(sigpipes[119], summers[1],0,29)]
                                   + [g.edge(sigpipes[120], summers[1],0,30)]
                                   + [g.edge(sigpipes[121], summers[1],0,31)]
                                   + [g.edge(sigpipes[122], summers[1],0,32)]
                                   + [g.edge(sigpipes[123], summers[1],0,33)]
                                   + [g.edge(sigpipes[124], summers[1],0,34)]
                                   + [g.edge(sigpipes[125], summers[1],0,35)]
                                   + [g.edge(sigpipes[126], summers[1],0,36)]
                                   + [g.edge(sigpipes[127], summers[1],0,37)]
                                   + [g.edge(sigpipes[128], summers[1],0,38)]
                                   + [g.edge(sigpipes[129], summers[1],0,39)]
                                   + [g.edge(sigpipes[130], summers[1],0,40)]
                                   + [g.edge(sigpipes[131], summers[1],0,41)]
                                   + [g.edge(sigpipes[132], summers[1],0,42)]
                                   + [g.edge(sigpipes[133], summers[1],0,43)]
                                   + [g.edge(sigpipes[134], summers[1],0,44)]
                                   + [g.edge(sigpipes[135], summers[1],0,45)]
                                   + [g.edge(sigpipes[136], summers[1],0,46)]
                                   + [g.edge(sigpipes[137], summers[1],0,47)]
                                   + [g.edge(sigpipes[138], summers[1],0,48)]
                                   + [g.edge(sigpipes[139], summers[1],0,49)]
                                   + [g.edge(sigpipes[140], summers[1],0,50)]
                                   + [g.edge(sigpipes[141], summers[1],0,51)]
                                   + [g.edge(sigpipes[142], summers[1],0,52)]
                                   + [g.edge(sigpipes[143], summers[1],0,53)]
                                   + [g.edge(sigpipes[144], summers[1],0,54)]
                                   + [g.edge(sigpipes[145], summers[1],0,55)]
                                   + [g.edge(sigpipes[146], summers[1],0,56)]
                                   + [g.edge(sigpipes[147], summers[1],0,57)]
                                   + [g.edge(sigpipes[148], summers[1],0,58)]
                                   + [g.edge(sigpipes[149], summers[1],0,59)]
                                   + [g.edge(sigpipes[150], summers[1],0,60)]
                                   + [g.edge(sigpipes[151], summers[1],0,61)]
                                   + [g.edge(sigpipes[152], summers[1],0,62)]
                                   + [g.edge(sigpipes[153], summers[1],0,63)]
                                   + [g.edge(sigpipes[154], summers[1],0,64)]
                                   + [g.edge(sigpipes[155], summers[1],0,65)]
                                   + [g.edge(sigpipes[156], summers[1],0,66)]
                                   + [g.edge(sigpipes[157], summers[1],0,67)]
                                   + [g.edge(sigpipes[158], summers[1],0,68)]
                                   + [g.edge(sigpipes[159], summers[1],0,69)]
                                   + [g.edge(sigpipes[160], summers[1],0,70)]
                                   + [g.edge(sigpipes[161], summers[1],0,71)]
                                   + [g.edge(sigpipes[162], summers[1],0,72)]
                                   + [g.edge(sigpipes[163], summers[1],0,73)]
                                   + [g.edge(sigpipes[164], summers[1],0,74)]
                                   + [g.edge(sigpipes[165], summers[1],0,75)]
                                   + [g.edge(sigpipes[166], summers[1],0,76)]
                                   + [g.edge(sigpipes[167], summers[1],0,77)]
                                   + [g.edge(sigpipes[168], summers[1],0,78)]
                                   + [g.edge(sigpipes[169], summers[1],0,79)]
                                   + [g.edge(sigpipes[170], summers[1],0,80)]
                                   + [g.edge(sigpipes[171], summers[1],0,81)]
                                   + [g.edge(sigpipes[172], summers[1],0,82)]
                                   + [g.edge(sigpipes[173], summers[1],0,83)]
                                   + [g.edge(sigpipes[174], summers[1],0,84)]
                                   + [g.edge(sigpipes[175], summers[1],0,85)]
                                   + [g.edge(sigpipes[176], summers[1],0,86)]
                                   + [g.edge(sigpipes[177], summers[1],0,87)]
                                   + [g.edge(sigpipes[178], summers[1],0,88)]
                                   + [g.edge(sigpipes[179], summers[1],0,89)]
                                   + [g.edge(sigpipes[180], summers[2],0,0)]
                                   + [g.edge(sigpipes[181], summers[2],0,1)]
                                   + [g.edge(sigpipes[182], summers[2],0,2)]
                                   + [g.edge(sigpipes[183], summers[2],0,3)]
                                   + [g.edge(sigpipes[184], summers[2],0,4)]
                                   + [g.edge(sigpipes[185], summers[2],0,5)]
                                   + [g.edge(sigpipes[186], summers[2],0,6)]
                                   + [g.edge(sigpipes[187], summers[2],0,7)]
                                   + [g.edge(sigpipes[188], summers[2],0,8)]
                                   + [g.edge(sigpipes[189], summers[2],0,9)]
                                   + [g.edge(sigpipes[190], summers[2],0,10)]
                                   + [g.edge(sigpipes[191], summers[2],0,11)]
                                   + [g.edge(sigpipes[192], summers[2],0,12)]
                                   + [g.edge(sigpipes[193], summers[2],0,13)]
                                   + [g.edge(sigpipes[194], summers[2],0,14)]
                                   + [g.edge(sigpipes[195], summers[2],0,15)]
                                   + [g.edge(sigpipes[196], summers[2],0,16)]
                                   + [g.edge(sigpipes[197], summers[2],0,17)]
                                   + [g.edge(sigpipes[198], summers[2],0,18)]
                                   + [g.edge(sigpipes[199], summers[2],0,19)]
                                   + [g.edge(sigpipes[200], summers[2],0,20)]
                                   + [g.edge(sigpipes[201], summers[2],0,21)]
                                   + [g.edge(sigpipes[202], summers[2],0,22)]
                                   + [g.edge(sigpipes[203], summers[2],0,23)]
                                   + [g.edge(sigpipes[204], summers[2],0,24)]
                                   + [g.edge(sigpipes[205], summers[2],0,25)]
                                   + [g.edge(sigpipes[206], summers[2],0,26)]
                                   + [g.edge(sigpipes[207], summers[2],0,27)]
                                   + [g.edge(sigpipes[208], summers[2],0,28)]
                                   + [g.edge(sigpipes[209], summers[2],0,29)]
                                   + [g.edge(sigpipes[210], summers[2],0,30)]
                                   + [g.edge(sigpipes[211], summers[2],0,31)]
                                   + [g.edge(sigpipes[212], summers[2],0,32)]
                                   + [g.edge(sigpipes[213], summers[2],0,33)]
                                   + [g.edge(sigpipes[214], summers[2],0,34)]
                                   + [g.edge(sigpipes[215], summers[2],0,35)]
                                   + [g.edge(sigpipes[216], summers[2],0,36)]
                                   + [g.edge(sigpipes[217], summers[2],0,37)]
                                   + [g.edge(sigpipes[218], summers[2],0,38)]
                                   + [g.edge(sigpipes[219], summers[2],0,39)]
                                   + [g.edge(sigpipes[220], summers[2],0,40)]
                                   + [g.edge(sigpipes[221], summers[2],0,41)]
                                   + [g.edge(sigpipes[222], summers[2],0,42)]
                                   + [g.edge(sigpipes[223], summers[2],0,43)]
                                   + [g.edge(sigpipes[224], summers[2],0,44)]
                                   + [g.edge(sigpipes[225], summers[2],0,45)]
                                   + [g.edge(sigpipes[226], summers[2],0,46)]
                                   + [g.edge(sigpipes[227], summers[2],0,47)]
                                   + [g.edge(sigpipes[228], summers[2],0,48)]
                                   + [g.edge(sigpipes[229], summers[2],0,49)]
                                   + [g.edge(sigpipes[230], summers[2],0,50)]
                                   + [g.edge(sigpipes[231], summers[2],0,51)]
                                   + [g.edge(sigpipes[232], summers[2],0,52)]
                                   + [g.edge(sigpipes[233], summers[2],0,53)]
                                   + [g.edge(sigpipes[234], summers[2],0,54)]
                                   + [g.edge(sigpipes[235], summers[2],0,55)]
                                   + [g.edge(sigpipes[236], summers[2],0,56)]
                                   + [g.edge(sigpipes[237], summers[2],0,57)]
                                   + [g.edge(sigpipes[238], summers[2],0,58)]
                                   + [g.edge(sigpipes[239], summers[2],0,59)]
                                   + [g.edge(sigpipes[240], summers[2],0,60)]
                                   + [g.edge(sigpipes[241], summers[2],0,61)]
                                   + [g.edge(sigpipes[242], summers[2],0,62)]
                                   + [g.edge(sigpipes[243], summers[2],0,63)]
                                   + [g.edge(sigpipes[244], summers[2],0,64)]
                                   + [g.edge(sigpipes[245], summers[2],0,65)]
                                   + [g.edge(sigpipes[246], summers[2],0,66)]
                                   + [g.edge(sigpipes[247], summers[2],0,67)]
                                   + [g.edge(sigpipes[248], summers[2],0,68)]
                                   + [g.edge(sigpipes[249], summers[2],0,69)]
                                   + [g.edge(sigpipes[250], summers[2],0,70)]
                                   + [g.edge(sigpipes[251], summers[2],0,71)]
                                   + [g.edge(sigpipes[252], summers[2],0,72)]
                                   + [g.edge(sigpipes[253], summers[2],0,73)]
                                   + [g.edge(sigpipes[254], summers[2],0,74)]
                                   + [g.edge(sigpipes[255], summers[2],0,75)]
                                   + [g.edge(sigpipes[256], summers[2],0,76)]
                                   + [g.edge(sigpipes[257], summers[2],0,77)]
                                   + [g.edge(sigpipes[258], summers[2],0,78)]
                                   + [g.edge(sigpipes[259], summers[2],0,79)]
                                   + [g.edge(sigpipes[260], summers[2],0,80)]
                                   + [g.edge(sigpipes[261], summers[2],0,81)]
                                   + [g.edge(sigpipes[262], summers[2],0,82)]
                                   + [g.edge(sigpipes[263], summers[2],0,83)]
                                   + [g.edge(sigpipes[264], summers[2],0,84)]
                                   + [g.edge(sigpipes[265], summers[2],0,85)]
                                   + [g.edge(sigpipes[266], summers[2],0,86)]
                                   + [g.edge(sigpipes[267], summers[2],0,87)]
                                   + [g.edge(sigpipes[268], summers[2],0,88)]
                                   + [g.edge(sigpipes[269], summers[2],0,89)]
                                   + [g.edge(sigpipes[270], summers[3],0,0)]
                                   + [g.edge(sigpipes[271], summers[3],0,1)]
                                   + [g.edge(sigpipes[272], summers[3],0,2)]
                                   + [g.edge(sigpipes[273], summers[3],0,3)]
                                   + [g.edge(sigpipes[274], summers[3],0,4)]
                                   + [g.edge(sigpipes[275], summers[3],0,5)]
                                   + [g.edge(sigpipes[276], summers[3],0,6)]
                                   + [g.edge(sigpipes[277], summers[3],0,7)]
                                   + [g.edge(sigpipes[278], summers[3],0,8)]
                                   + [g.edge(sigpipes[279], summers[3],0,9)]
                                   + [g.edge(sigpipes[280], summers[3],0,10)]
                                   + [g.edge(sigpipes[281], summers[3],0,11)]
                                   + [g.edge(sigpipes[282], summers[3],0,12)]
                                   + [g.edge(sigpipes[283], summers[3],0,13)]
                                   + [g.edge(sigpipes[284], summers[3],0,14)]
                                   + [g.edge(sigpipes[285], summers[3],0,15)]
                                   + [g.edge(sigpipes[286], summers[3],0,16)]
                                   + [g.edge(sigpipes[287], summers[3],0,17)]
                                   + [g.edge(sigpipes[288], summers[3],0,18)]
                                   + [g.edge(sigpipes[289], summers[3],0,19)]
                                   + [g.edge(sigpipes[290], summers[3],0,20)]
                                   + [g.edge(sigpipes[291], summers[3],0,21)]
                                   + [g.edge(sigpipes[292], summers[3],0,22)]
                                   + [g.edge(sigpipes[293], summers[3],0,23)]
                                   + [g.edge(sigpipes[294], summers[3],0,24)]
                                   + [g.edge(sigpipes[295], summers[3],0,25)]
                                   + [g.edge(sigpipes[296], summers[3],0,26)]
                                   + [g.edge(sigpipes[297], summers[3],0,27)]
                                   + [g.edge(sigpipes[298], summers[3],0,28)]
                                   + [g.edge(sigpipes[299], summers[3],0,29)]
                                   + [g.edge(sigpipes[300], summers[3],0,30)]
                                   + [g.edge(sigpipes[301], summers[3],0,31)]
                                   + [g.edge(sigpipes[302], summers[3],0,32)]
                                   + [g.edge(sigpipes[303], summers[3],0,33)]
                                   + [g.edge(sigpipes[304], summers[3],0,34)]
                                   + [g.edge(sigpipes[305], summers[3],0,35)]
                                   + [g.edge(sigpipes[306], summers[3],0,36)]
                                   + [g.edge(sigpipes[307], summers[3],0,37)]
                                   + [g.edge(sigpipes[308], summers[3],0,38)]
                                   + [g.edge(sigpipes[309], summers[3],0,39)]
                                   + [g.edge(sigpipes[310], summers[3],0,40)]
                                   + [g.edge(sigpipes[311], summers[3],0,41)]
                                   + [g.edge(sigpipes[312], summers[3],0,42)]
                                   + [g.edge(sigpipes[313], summers[3],0,43)]
                                   + [g.edge(sigpipes[314], summers[3],0,44)]
                                   + [g.edge(sigpipes[315], summers[3],0,45)]
                                   + [g.edge(sigpipes[316], summers[3],0,46)]
                                   + [g.edge(sigpipes[317], summers[3],0,47)]
                                   + [g.edge(sigpipes[318], summers[3],0,48)]
                                   + [g.edge(sigpipes[319], summers[3],0,49)]
                                   + [g.edge(sigpipes[320], summers[3],0,50)]
                                   + [g.edge(sigpipes[321], summers[3],0,51)]
                                   + [g.edge(sigpipes[322], summers[3],0,52)]
                                   + [g.edge(sigpipes[323], summers[3],0,53)]
                                   + [g.edge(sigpipes[324], summers[3],0,54)]
                                   + [g.edge(sigpipes[325], summers[3],0,55)]
                                   + [g.edge(sigpipes[326], summers[3],0,56)]
                                   + [g.edge(sigpipes[327], summers[3],0,57)]
                                   + [g.edge(sigpipes[328], summers[3],0,58)]
                                   + [g.edge(sigpipes[329], summers[3],0,59)]
                                   + [g.edge(sigpipes[330], summers[3],0,60)]
                                   + [g.edge(sigpipes[331], summers[3],0,61)]
                                   + [g.edge(sigpipes[332], summers[3],0,62)]
                                   + [g.edge(sigpipes[333], summers[3],0,63)]
                                   + [g.edge(sigpipes[334], summers[3],0,64)]
                                   + [g.edge(sigpipes[335], summers[3],0,65)]
                                   + [g.edge(sigpipes[336], summers[3],0,66)]
                                   + [g.edge(sigpipes[337], summers[3],0,67)]
                                   + [g.edge(sigpipes[338], summers[3],0,68)]
                                   + [g.edge(sigpipes[339], summers[3],0,69)]
                                   + [g.edge(sigpipes[340], summers[3],0,70)]
                                   + [g.edge(sigpipes[341], summers[3],0,71)]
                                   + [g.edge(sigpipes[342], summers[3],0,72)]
                                   + [g.edge(sigpipes[343], summers[3],0,73)]
                                   + [g.edge(sigpipes[344], summers[3],0,74)]
                                   + [g.edge(sigpipes[345], summers[3],0,75)]
                                   + [g.edge(sigpipes[346], summers[3],0,76)]
                                   + [g.edge(sigpipes[347], summers[3],0,77)]
                                   + [g.edge(sigpipes[348], summers[3],0,78)]
                                   + [g.edge(sigpipes[349], summers[3],0,79)]
                                   + [g.edge(sigpipes[350], summers[3],0,80)]
                                   + [g.edge(sigpipes[351], summers[3],0,81)]
                                   + [g.edge(sigpipes[352], summers[3],0,82)]
                                   + [g.edge(sigpipes[353], summers[3],0,83)]
                                   + [g.edge(sigpipes[354], summers[3],0,84)]
                                   + [g.edge(sigpipes[355], summers[3],0,85)]
                                   + [g.edge(sigpipes[356], summers[3],0,86)]
                                   + [g.edge(sigpipes[357], summers[3],0,87)]
                                   + [g.edge(sigpipes[358], summers[3],0,88)]
                                   + [g.edge(sigpipes[359], summers[3],0,89)]
                                 // connecting summer and the operator pipelines
                                 + [g.edge(summers[n], actpipes[n]) for n in std.range(0,faninmult-1)],
                                 name=name),

        local signal = g.intern(innodes=sigpipes,
                                 outnodes=actpipes,
                                 centernodes=summers,
                                 edges= 
                                 // connecting signal and summer
                                   + [g.edge(sigpipes[0], summers[0],0,0)]
                                   + [g.edge(sigpipes[1], summers[0],0,1)]
                                   + [g.edge(sigpipes[2], summers[0],0,2)]
                                   + [g.edge(sigpipes[3], summers[0],0,3)]
                                   + [g.edge(sigpipes[4], summers[0],0,4)]
                                   + [g.edge(sigpipes[5], summers[0],0,5)]
                                   + [g.edge(sigpipes[6], summers[0],0,6)]
                                   + [g.edge(sigpipes[7], summers[0],0,7)]
                                   + [g.edge(sigpipes[8], summers[0],0,8)]
                                   + [g.edge(sigpipes[9], summers[0],0,9)]
                                   + [g.edge(sigpipes[10], summers[0],0,10)]
                                   + [g.edge(sigpipes[11], summers[0],0,11)]
                                   + [g.edge(sigpipes[12], summers[0],0,12)]
                                   + [g.edge(sigpipes[13], summers[0],0,13)]
                                   + [g.edge(sigpipes[14], summers[0],0,14)]
                                   + [g.edge(sigpipes[15], summers[0],0,15)]
                                   + [g.edge(sigpipes[16], summers[0],0,16)]
                                   + [g.edge(sigpipes[17], summers[0],0,17)]
                                   + [g.edge(sigpipes[18], summers[0],0,18)]
                                   + [g.edge(sigpipes[19], summers[0],0,19)]
                                   + [g.edge(sigpipes[20], summers[0],0,20)]
                                   + [g.edge(sigpipes[21], summers[0],0,21)]
                                   + [g.edge(sigpipes[22], summers[0],0,22)]
                                   + [g.edge(sigpipes[23], summers[0],0,23)]
                                   + [g.edge(sigpipes[24], summers[0],0,24)]
                                   + [g.edge(sigpipes[25], summers[0],0,25)]
                                   + [g.edge(sigpipes[26], summers[0],0,26)]
                                   + [g.edge(sigpipes[27], summers[0],0,27)]
                                   + [g.edge(sigpipes[28], summers[0],0,28)]
                                   + [g.edge(sigpipes[29], summers[0],0,29)]
                                   + [g.edge(sigpipes[30], summers[0],0,30)]
                                   + [g.edge(sigpipes[31], summers[0],0,31)]
                                   + [g.edge(sigpipes[32], summers[0],0,32)]
                                   + [g.edge(sigpipes[33], summers[0],0,33)]
                                   + [g.edge(sigpipes[34], summers[0],0,34)]
                                   + [g.edge(sigpipes[35], summers[0],0,35)]
                                   + [g.edge(sigpipes[36], summers[0],0,36)]
                                   + [g.edge(sigpipes[37], summers[0],0,37)]
                                   + [g.edge(sigpipes[38], summers[0],0,38)]
                                   + [g.edge(sigpipes[39], summers[0],0,39)]
                                   + [g.edge(sigpipes[40], summers[0],0,40)]
                                   + [g.edge(sigpipes[41], summers[0],0,41)]
                                   + [g.edge(sigpipes[42], summers[0],0,42)]
                                   + [g.edge(sigpipes[43], summers[0],0,43)]
                                   + [g.edge(sigpipes[44], summers[0],0,44)]
                                   + [g.edge(sigpipes[45], summers[0],0,45)]
                                   + [g.edge(sigpipes[46], summers[0],0,46)]
                                   + [g.edge(sigpipes[47], summers[0],0,47)]
                                   + [g.edge(sigpipes[48], summers[0],0,48)]
                                   + [g.edge(sigpipes[49], summers[0],0,49)]
                                   + [g.edge(sigpipes[50], summers[0],0,50)]
                                   + [g.edge(sigpipes[51], summers[0],0,51)]
                                   + [g.edge(sigpipes[52], summers[0],0,52)]
                                   + [g.edge(sigpipes[53], summers[0],0,53)]
                                   + [g.edge(sigpipes[54], summers[0],0,54)]
                                   + [g.edge(sigpipes[55], summers[0],0,55)]
                                   + [g.edge(sigpipes[56], summers[0],0,56)]
                                   + [g.edge(sigpipes[57], summers[0],0,57)]
                                   + [g.edge(sigpipes[58], summers[0],0,58)]
                                   + [g.edge(sigpipes[59], summers[0],0,59)]
                                   + [g.edge(sigpipes[60], summers[0],0,60)]
                                   + [g.edge(sigpipes[61], summers[0],0,61)]
                                   + [g.edge(sigpipes[62], summers[0],0,62)]
                                   + [g.edge(sigpipes[63], summers[0],0,63)]
                                   + [g.edge(sigpipes[64], summers[0],0,64)]
                                   + [g.edge(sigpipes[65], summers[0],0,65)]
                                   + [g.edge(sigpipes[66], summers[0],0,66)]
                                   + [g.edge(sigpipes[67], summers[0],0,67)]
                                   + [g.edge(sigpipes[68], summers[0],0,68)]
                                   + [g.edge(sigpipes[69], summers[0],0,69)]
                                   + [g.edge(sigpipes[70], summers[0],0,70)]
                                   + [g.edge(sigpipes[71], summers[0],0,71)]
                                   + [g.edge(sigpipes[72], summers[0],0,72)]
                                   + [g.edge(sigpipes[73], summers[0],0,73)]
                                   + [g.edge(sigpipes[74], summers[0],0,74)]
                                   + [g.edge(sigpipes[75], summers[0],0,75)]
                                   + [g.edge(sigpipes[76], summers[0],0,76)]
                                   + [g.edge(sigpipes[77], summers[0],0,77)]
                                   + [g.edge(sigpipes[78], summers[0],0,78)]
                                   + [g.edge(sigpipes[79], summers[0],0,79)]
                                   + [g.edge(sigpipes[80], summers[0],0,80)]
                                   + [g.edge(sigpipes[81], summers[0],0,81)]
                                   + [g.edge(sigpipes[82], summers[0],0,82)]
                                   + [g.edge(sigpipes[83], summers[0],0,83)]
                                   + [g.edge(sigpipes[84], summers[0],0,84)]
                                   + [g.edge(sigpipes[85], summers[0],0,85)]
                                   + [g.edge(sigpipes[86], summers[0],0,86)]
                                   + [g.edge(sigpipes[87], summers[0],0,87)]
                                   + [g.edge(sigpipes[88], summers[0],0,88)]
                                   + [g.edge(sigpipes[89], summers[0],0,89)]
                                   + [g.edge(sigpipes[90], summers[1],0,0)]
                                   + [g.edge(sigpipes[91], summers[1],0,1)]
                                   + [g.edge(sigpipes[92], summers[1],0,2)]
                                   + [g.edge(sigpipes[93], summers[1],0,3)]
                                   + [g.edge(sigpipes[94], summers[1],0,4)]
                                   + [g.edge(sigpipes[95], summers[1],0,5)]
                                   + [g.edge(sigpipes[96], summers[1],0,6)]
                                   + [g.edge(sigpipes[97], summers[1],0,7)]
                                   + [g.edge(sigpipes[98], summers[1],0,8)]
                                   + [g.edge(sigpipes[99], summers[1],0,9)]
                                   + [g.edge(sigpipes[100], summers[1],0,10)]
                                   + [g.edge(sigpipes[101], summers[1],0,11)]
                                   + [g.edge(sigpipes[102], summers[1],0,12)]
                                   + [g.edge(sigpipes[103], summers[1],0,13)]
                                   + [g.edge(sigpipes[104], summers[1],0,14)]
                                   + [g.edge(sigpipes[105], summers[1],0,15)]
                                   + [g.edge(sigpipes[106], summers[1],0,16)]
                                   + [g.edge(sigpipes[107], summers[1],0,17)]
                                   + [g.edge(sigpipes[108], summers[1],0,18)]
                                   + [g.edge(sigpipes[109], summers[1],0,19)]
                                   + [g.edge(sigpipes[110], summers[1],0,20)]
                                   + [g.edge(sigpipes[111], summers[1],0,21)]
                                   + [g.edge(sigpipes[112], summers[1],0,22)]
                                   + [g.edge(sigpipes[113], summers[1],0,23)]
                                   + [g.edge(sigpipes[114], summers[1],0,24)]
                                   + [g.edge(sigpipes[115], summers[1],0,25)]
                                   + [g.edge(sigpipes[116], summers[1],0,26)]
                                   + [g.edge(sigpipes[117], summers[1],0,27)]
                                   + [g.edge(sigpipes[118], summers[1],0,28)]
                                   + [g.edge(sigpipes[119], summers[1],0,29)]
                                   + [g.edge(sigpipes[120], summers[1],0,30)]
                                   + [g.edge(sigpipes[121], summers[1],0,31)]
                                   + [g.edge(sigpipes[122], summers[1],0,32)]
                                   + [g.edge(sigpipes[123], summers[1],0,33)]
                                   + [g.edge(sigpipes[124], summers[1],0,34)]
                                   + [g.edge(sigpipes[125], summers[1],0,35)]
                                   + [g.edge(sigpipes[126], summers[1],0,36)]
                                   + [g.edge(sigpipes[127], summers[1],0,37)]
                                   + [g.edge(sigpipes[128], summers[1],0,38)]
                                   + [g.edge(sigpipes[129], summers[1],0,39)]
                                   + [g.edge(sigpipes[130], summers[1],0,40)]
                                   + [g.edge(sigpipes[131], summers[1],0,41)]
                                   + [g.edge(sigpipes[132], summers[1],0,42)]
                                   + [g.edge(sigpipes[133], summers[1],0,43)]
                                   + [g.edge(sigpipes[134], summers[1],0,44)]
                                   + [g.edge(sigpipes[135], summers[1],0,45)]
                                   + [g.edge(sigpipes[136], summers[1],0,46)]
                                   + [g.edge(sigpipes[137], summers[1],0,47)]
                                   + [g.edge(sigpipes[138], summers[1],0,48)]
                                   + [g.edge(sigpipes[139], summers[1],0,49)]
                                   + [g.edge(sigpipes[140], summers[1],0,50)]
                                   + [g.edge(sigpipes[141], summers[1],0,51)]
                                   + [g.edge(sigpipes[142], summers[1],0,52)]
                                   + [g.edge(sigpipes[143], summers[1],0,53)]
                                   + [g.edge(sigpipes[144], summers[1],0,54)]
                                   + [g.edge(sigpipes[145], summers[1],0,55)]
                                   + [g.edge(sigpipes[146], summers[1],0,56)]
                                   + [g.edge(sigpipes[147], summers[1],0,57)]
                                   + [g.edge(sigpipes[148], summers[1],0,58)]
                                   + [g.edge(sigpipes[149], summers[1],0,59)]
                                   + [g.edge(sigpipes[150], summers[1],0,60)]
                                   + [g.edge(sigpipes[151], summers[1],0,61)]
                                   + [g.edge(sigpipes[152], summers[1],0,62)]
                                   + [g.edge(sigpipes[153], summers[1],0,63)]
                                   + [g.edge(sigpipes[154], summers[1],0,64)]
                                   + [g.edge(sigpipes[155], summers[1],0,65)]
                                   + [g.edge(sigpipes[156], summers[1],0,66)]
                                   + [g.edge(sigpipes[157], summers[1],0,67)]
                                   + [g.edge(sigpipes[158], summers[1],0,68)]
                                   + [g.edge(sigpipes[159], summers[1],0,69)]
                                   + [g.edge(sigpipes[160], summers[1],0,70)]
                                   + [g.edge(sigpipes[161], summers[1],0,71)]
                                   + [g.edge(sigpipes[162], summers[1],0,72)]
                                   + [g.edge(sigpipes[163], summers[1],0,73)]
                                   + [g.edge(sigpipes[164], summers[1],0,74)]
                                   + [g.edge(sigpipes[165], summers[1],0,75)]
                                   + [g.edge(sigpipes[166], summers[1],0,76)]
                                   + [g.edge(sigpipes[167], summers[1],0,77)]
                                   + [g.edge(sigpipes[168], summers[1],0,78)]
                                   + [g.edge(sigpipes[169], summers[1],0,79)]
                                   + [g.edge(sigpipes[170], summers[1],0,80)]
                                   + [g.edge(sigpipes[171], summers[1],0,81)]
                                   + [g.edge(sigpipes[172], summers[1],0,82)]
                                   + [g.edge(sigpipes[173], summers[1],0,83)]
                                   + [g.edge(sigpipes[174], summers[1],0,84)]
                                   + [g.edge(sigpipes[175], summers[1],0,85)]
                                   + [g.edge(sigpipes[176], summers[1],0,86)]
                                   + [g.edge(sigpipes[177], summers[1],0,87)]
                                   + [g.edge(sigpipes[178], summers[1],0,88)]
                                   + [g.edge(sigpipes[179], summers[1],0,89)]
                                   + [g.edge(sigpipes[180], summers[2],0,0)]
                                   + [g.edge(sigpipes[181], summers[2],0,1)]
                                   + [g.edge(sigpipes[182], summers[2],0,2)]
                                   + [g.edge(sigpipes[183], summers[2],0,3)]
                                   + [g.edge(sigpipes[184], summers[2],0,4)]
                                   + [g.edge(sigpipes[185], summers[2],0,5)]
                                   + [g.edge(sigpipes[186], summers[2],0,6)]
                                   + [g.edge(sigpipes[187], summers[2],0,7)]
                                   + [g.edge(sigpipes[188], summers[2],0,8)]
                                   + [g.edge(sigpipes[189], summers[2],0,9)]
                                   + [g.edge(sigpipes[190], summers[2],0,10)]
                                   + [g.edge(sigpipes[191], summers[2],0,11)]
                                   + [g.edge(sigpipes[192], summers[2],0,12)]
                                   + [g.edge(sigpipes[193], summers[2],0,13)]
                                   + [g.edge(sigpipes[194], summers[2],0,14)]
                                   + [g.edge(sigpipes[195], summers[2],0,15)]
                                   + [g.edge(sigpipes[196], summers[2],0,16)]
                                   + [g.edge(sigpipes[197], summers[2],0,17)]
                                   + [g.edge(sigpipes[198], summers[2],0,18)]
                                   + [g.edge(sigpipes[199], summers[2],0,19)]
                                   + [g.edge(sigpipes[200], summers[2],0,20)]
                                   + [g.edge(sigpipes[201], summers[2],0,21)]
                                   + [g.edge(sigpipes[202], summers[2],0,22)]
                                   + [g.edge(sigpipes[203], summers[2],0,23)]
                                   + [g.edge(sigpipes[204], summers[2],0,24)]
                                   + [g.edge(sigpipes[205], summers[2],0,25)]
                                   + [g.edge(sigpipes[206], summers[2],0,26)]
                                   + [g.edge(sigpipes[207], summers[2],0,27)]
                                   + [g.edge(sigpipes[208], summers[2],0,28)]
                                   + [g.edge(sigpipes[209], summers[2],0,29)]
                                   + [g.edge(sigpipes[210], summers[2],0,30)]
                                   + [g.edge(sigpipes[211], summers[2],0,31)]
                                   + [g.edge(sigpipes[212], summers[2],0,32)]
                                   + [g.edge(sigpipes[213], summers[2],0,33)]
                                   + [g.edge(sigpipes[214], summers[2],0,34)]
                                   + [g.edge(sigpipes[215], summers[2],0,35)]
                                   + [g.edge(sigpipes[216], summers[2],0,36)]
                                   + [g.edge(sigpipes[217], summers[2],0,37)]
                                   + [g.edge(sigpipes[218], summers[2],0,38)]
                                   + [g.edge(sigpipes[219], summers[2],0,39)]
                                   + [g.edge(sigpipes[220], summers[2],0,40)]
                                   + [g.edge(sigpipes[221], summers[2],0,41)]
                                   + [g.edge(sigpipes[222], summers[2],0,42)]
                                   + [g.edge(sigpipes[223], summers[2],0,43)]
                                   + [g.edge(sigpipes[224], summers[2],0,44)]
                                   + [g.edge(sigpipes[225], summers[2],0,45)]
                                   + [g.edge(sigpipes[226], summers[2],0,46)]
                                   + [g.edge(sigpipes[227], summers[2],0,47)]
                                   + [g.edge(sigpipes[228], summers[2],0,48)]
                                   + [g.edge(sigpipes[229], summers[2],0,49)]
                                   + [g.edge(sigpipes[230], summers[2],0,50)]
                                   + [g.edge(sigpipes[231], summers[2],0,51)]
                                   + [g.edge(sigpipes[232], summers[2],0,52)]
                                   + [g.edge(sigpipes[233], summers[2],0,53)]
                                   + [g.edge(sigpipes[234], summers[2],0,54)]
                                   + [g.edge(sigpipes[235], summers[2],0,55)]
                                   + [g.edge(sigpipes[236], summers[2],0,56)]
                                   + [g.edge(sigpipes[237], summers[2],0,57)]
                                   + [g.edge(sigpipes[238], summers[2],0,58)]
                                   + [g.edge(sigpipes[239], summers[2],0,59)]
                                   + [g.edge(sigpipes[240], summers[2],0,60)]
                                   + [g.edge(sigpipes[241], summers[2],0,61)]
                                   + [g.edge(sigpipes[242], summers[2],0,62)]
                                   + [g.edge(sigpipes[243], summers[2],0,63)]
                                   + [g.edge(sigpipes[244], summers[2],0,64)]
                                   + [g.edge(sigpipes[245], summers[2],0,65)]
                                   + [g.edge(sigpipes[246], summers[2],0,66)]
                                   + [g.edge(sigpipes[247], summers[2],0,67)]
                                   + [g.edge(sigpipes[248], summers[2],0,68)]
                                   + [g.edge(sigpipes[249], summers[2],0,69)]
                                   + [g.edge(sigpipes[250], summers[2],0,70)]
                                   + [g.edge(sigpipes[251], summers[2],0,71)]
                                   + [g.edge(sigpipes[252], summers[2],0,72)]
                                   + [g.edge(sigpipes[253], summers[2],0,73)]
                                   + [g.edge(sigpipes[254], summers[2],0,74)]
                                   + [g.edge(sigpipes[255], summers[2],0,75)]
                                   + [g.edge(sigpipes[256], summers[2],0,76)]
                                   + [g.edge(sigpipes[257], summers[2],0,77)]
                                   + [g.edge(sigpipes[258], summers[2],0,78)]
                                   + [g.edge(sigpipes[259], summers[2],0,79)]
                                   + [g.edge(sigpipes[260], summers[2],0,80)]
                                   + [g.edge(sigpipes[261], summers[2],0,81)]
                                   + [g.edge(sigpipes[262], summers[2],0,82)]
                                   + [g.edge(sigpipes[263], summers[2],0,83)]
                                   + [g.edge(sigpipes[264], summers[2],0,84)]
                                   + [g.edge(sigpipes[265], summers[2],0,85)]
                                   + [g.edge(sigpipes[266], summers[2],0,86)]
                                   + [g.edge(sigpipes[267], summers[2],0,87)]
                                   + [g.edge(sigpipes[268], summers[2],0,88)]
                                   + [g.edge(sigpipes[269], summers[2],0,89)]
                                   + [g.edge(sigpipes[270], summers[3],0,0)]
                                   + [g.edge(sigpipes[271], summers[3],0,1)]
                                   + [g.edge(sigpipes[272], summers[3],0,2)]
                                   + [g.edge(sigpipes[273], summers[3],0,3)]
                                   + [g.edge(sigpipes[274], summers[3],0,4)]
                                   + [g.edge(sigpipes[275], summers[3],0,5)]
                                   + [g.edge(sigpipes[276], summers[3],0,6)]
                                   + [g.edge(sigpipes[277], summers[3],0,7)]
                                   + [g.edge(sigpipes[278], summers[3],0,8)]
                                   + [g.edge(sigpipes[279], summers[3],0,9)]
                                   + [g.edge(sigpipes[280], summers[3],0,10)]
                                   + [g.edge(sigpipes[281], summers[3],0,11)]
                                   + [g.edge(sigpipes[282], summers[3],0,12)]
                                   + [g.edge(sigpipes[283], summers[3],0,13)]
                                   + [g.edge(sigpipes[284], summers[3],0,14)]
                                   + [g.edge(sigpipes[285], summers[3],0,15)]
                                   + [g.edge(sigpipes[286], summers[3],0,16)]
                                   + [g.edge(sigpipes[287], summers[3],0,17)]
                                   + [g.edge(sigpipes[288], summers[3],0,18)]
                                   + [g.edge(sigpipes[289], summers[3],0,19)]
                                   + [g.edge(sigpipes[290], summers[3],0,20)]
                                   + [g.edge(sigpipes[291], summers[3],0,21)]
                                   + [g.edge(sigpipes[292], summers[3],0,22)]
                                   + [g.edge(sigpipes[293], summers[3],0,23)]
                                   + [g.edge(sigpipes[294], summers[3],0,24)]
                                   + [g.edge(sigpipes[295], summers[3],0,25)]
                                   + [g.edge(sigpipes[296], summers[3],0,26)]
                                   + [g.edge(sigpipes[297], summers[3],0,27)]
                                   + [g.edge(sigpipes[298], summers[3],0,28)]
                                   + [g.edge(sigpipes[299], summers[3],0,29)]
                                   + [g.edge(sigpipes[300], summers[3],0,30)]
                                   + [g.edge(sigpipes[301], summers[3],0,31)]
                                   + [g.edge(sigpipes[302], summers[3],0,32)]
                                   + [g.edge(sigpipes[303], summers[3],0,33)]
                                   + [g.edge(sigpipes[304], summers[3],0,34)]
                                   + [g.edge(sigpipes[305], summers[3],0,35)]
                                   + [g.edge(sigpipes[306], summers[3],0,36)]
                                   + [g.edge(sigpipes[307], summers[3],0,37)]
                                   + [g.edge(sigpipes[308], summers[3],0,38)]
                                   + [g.edge(sigpipes[309], summers[3],0,39)]
                                   + [g.edge(sigpipes[310], summers[3],0,40)]
                                   + [g.edge(sigpipes[311], summers[3],0,41)]
                                   + [g.edge(sigpipes[312], summers[3],0,42)]
                                   + [g.edge(sigpipes[313], summers[3],0,43)]
                                   + [g.edge(sigpipes[314], summers[3],0,44)]
                                   + [g.edge(sigpipes[315], summers[3],0,45)]
                                   + [g.edge(sigpipes[316], summers[3],0,46)]
                                   + [g.edge(sigpipes[317], summers[3],0,47)]
                                   + [g.edge(sigpipes[318], summers[3],0,48)]
                                   + [g.edge(sigpipes[319], summers[3],0,49)]
                                   + [g.edge(sigpipes[320], summers[3],0,50)]
                                   + [g.edge(sigpipes[321], summers[3],0,51)]
                                   + [g.edge(sigpipes[322], summers[3],0,52)]
                                   + [g.edge(sigpipes[323], summers[3],0,53)]
                                   + [g.edge(sigpipes[324], summers[3],0,54)]
                                   + [g.edge(sigpipes[325], summers[3],0,55)]
                                   + [g.edge(sigpipes[326], summers[3],0,56)]
                                   + [g.edge(sigpipes[327], summers[3],0,57)]
                                   + [g.edge(sigpipes[328], summers[3],0,58)]
                                   + [g.edge(sigpipes[329], summers[3],0,59)]
                                   + [g.edge(sigpipes[330], summers[3],0,60)]
                                   + [g.edge(sigpipes[331], summers[3],0,61)]
                                   + [g.edge(sigpipes[332], summers[3],0,62)]
                                   + [g.edge(sigpipes[333], summers[3],0,63)]
                                   + [g.edge(sigpipes[334], summers[3],0,64)]
                                   + [g.edge(sigpipes[335], summers[3],0,65)]
                                   + [g.edge(sigpipes[336], summers[3],0,66)]
                                   + [g.edge(sigpipes[337], summers[3],0,67)]
                                   + [g.edge(sigpipes[338], summers[3],0,68)]
                                   + [g.edge(sigpipes[339], summers[3],0,69)]
                                   + [g.edge(sigpipes[340], summers[3],0,70)]
                                   + [g.edge(sigpipes[341], summers[3],0,71)]
                                   + [g.edge(sigpipes[342], summers[3],0,72)]
                                   + [g.edge(sigpipes[343], summers[3],0,73)]
                                   + [g.edge(sigpipes[344], summers[3],0,74)]
                                   + [g.edge(sigpipes[345], summers[3],0,75)]
                                   + [g.edge(sigpipes[346], summers[3],0,76)]
                                   + [g.edge(sigpipes[347], summers[3],0,77)]
                                   + [g.edge(sigpipes[348], summers[3],0,78)]
                                   + [g.edge(sigpipes[349], summers[3],0,79)]
                                   + [g.edge(sigpipes[350], summers[3],0,80)]
                                   + [g.edge(sigpipes[351], summers[3],0,81)]
                                   + [g.edge(sigpipes[352], summers[3],0,82)]
                                   + [g.edge(sigpipes[353], summers[3],0,83)]
                                   + [g.edge(sigpipes[354], summers[3],0,84)]
                                   + [g.edge(sigpipes[355], summers[3],0,85)]
                                   + [g.edge(sigpipes[356], summers[3],0,86)]
                                   + [g.edge(sigpipes[357], summers[3],0,87)]
                                   + [g.edge(sigpipes[358], summers[3],0,88)]
                                   + [g.edge(sigpipes[359], summers[3],0,89)]
                                 // connecting summer and the operator pipelines
                                 + [g.edge(summers[n], actpipes[n]) for n in std.range(0,faninmult-1)],
                                 name=name),


        ret: g.intern(innodes=[fanout],
                      outnodes=[fanin],
                      centernodes=[drift],
                      edges=
                      [g.edge(fanout, driftpipes[n], n, 0) for n in std.range(0, fanoutmult-1)] +
                      [g.edge(drift, fanin, n, n) for n in std.range(0, faninmult-1)],
                      name=name),
    }.ret,





  // Build a fanout-[pipelines]-fanin graph.  pipelines is a list of
  // pnode objects, one for each spine of the fan.
  fanpipe:: function(fout, pipelines, fin, name='fanpipe', outtags=[], fout_tag_rules=[], fin_tag_rules=[]) {

    local fanmult = std.length(pipelines),
    local fannums = std.range(0, fanmult - 1),

    local fanout = g.pnode({
      type: fout,
      name: name,
      data: {
        multiplicity: fanmult,
        tag_rules: fout_tag_rules,
      },
    }, nin=1, nout=fanmult),


    local fanin = g.pnode({
      type: fin,
      name: name,
      data: {
        multiplicity: fanmult,
        tag_rules: fin_tag_rules,
        tags: outtags,
      },
    }, nin=fanmult, nout=1),

    ret: g.intern(innodes=[fanout],
                  outnodes=[fanin],
                  centernodes=pipelines,
                  edges=
                  [g.edge(fanout, pipelines[n], n, 0) for n in std.range(0, fanmult - 1)] +
                  [g.edge(pipelines[n], fanin, 0, n) for n in std.range(0, fanmult - 1)],
                  name=name),
  }.ret,

      multifanpipe :: function( fout, pipelines, fin,
    fout_nnodes=[1,8,16], fout_multi=[8,2,7],
    fin_nnodes=[1,8,16], fin_multi=[8,2,7],
    name='multifanpipe', outtags=[], tag_rules=null ) {
        local fout_nlayers = std.length(fout_multi),
        assert fout_nlayers >= 2 : "fout_nlayers should be >= 2",
        local fin_nlayers = std.length(fin_multi),
        assert fin_nlayers >= 2 : "fin_nlayers should be >= 2",
        local npipe = std.length(pipelines),
        assert npipe == fout_nnodes[fout_nlayers-1]*fout_multi[fout_nlayers-1] :
            "fout layout error npipe=%d, "%npipe + "fout=%d"%(fout_nnodes[fout_nlayers-1]*fout_multi[fout_nlayers-1]),
        assert npipe == fin_nnodes[fin_nlayers-1]*fin_multi[fin_nlayers-1] :
            "fin layout error npipe=%d, "%npipe + "fin=%d"%(fin_nnodes[fin_nlayers-1]*fin_multi[fin_nlayers-1]),

        // function to create nodes for one layer
        local fout_layer(ilayer,nnodes,nmulti) = {
            ret : [
                g.pnode({
                    type: fout,
                    name: name+"_fout_%d"%ilayer + "_%d"%inode,
                    data: {
                        multiplicity: nmulti,
                        tag_rules: [],
                    }}, nin=1, nout=nmulti
                ) for inode in std.range(0,nnodes-1)],
        }.ret,
        // nodes for all layers
        local fout_layers = [
            fout_layer(ilayer,
            fout_nnodes[ilayer],
            fout_multi[ilayer])
            for ilayer in std.range(0,fout_nlayers-1)
        ],
        // make edges to make a combo node
        local fout_node = g.intern(
            innodes = fout_layers[0],
            centernodes = if fout_nlayers == 2 then [] else std.flattenArrays([fout_layers[i] for i in std.range(1,fout_nlayers-2)]),
            outnodes = fout_layers[fout_nlayers-1],
            edges = std.flattenArrays(
                [
                    [
                        g.edge(
                            fout_layers[ilayer-1][std.floor(inode/fout_multi[ilayer-1])],
                            fout_layers[ilayer][inode],
                            inode%fout_multi[ilayer-1],
                            0)
                    for inode in std.range(0,fout_nnodes[ilayer]-1)]
                for ilayer in std.range(1,fout_nlayers-1)])
        ),

        // similarly build the multi-layer fan in combo node
        // note the backward layer counting
        local fin_layer(ilayer,nnodes,nmulti) = {
            ret : [
                g.pnode({
                    type: fin,
                    name: name+"_fin_%d"%ilayer + "_%d"%inode,
                    data: {
                        multiplicity: nmulti,
                        tags: outtags,
                        tag_rules: [tag_rules for irule in std.range(0,nmulti-1)],
                    }}, nin=nmulti, nout=1
                ) for inode in std.range(0,nnodes-1)],
        }.ret,
        local fin_layers = [
            fin_layer(ilayer,
            fin_nnodes[ilayer],
            fin_multi[ilayer])
            for ilayer in std.range(0,fin_nlayers-1)
        ],
        local fin_node = g.intern(
            innodes = fin_layers[fin_nlayers-1],
            centernodes = if fin_nlayers == 2 then [] else std.flattenArrays([fin_layers[i] for i in std.range(1,fout_nlayers-2)]),
            outnodes = fin_layers[0],
            edges = std.flattenArrays(
                [
                    [
                        g.edge(
                            fin_layers[ilayer][inode],
                            fin_layers[ilayer-1][std.floor(inode/fin_multi[ilayer-1])],
                            0,
                            inode%fin_multi[ilayer-1])
                    for inode in std.range(0,fin_nnodes[ilayer]-1)]
                for ilayer in std.range(1,fin_nlayers-1)])
        ),

        // connect comb_fan_out-piples-combo_fan_in
        ret : g.intern(
            innodes = [fout_node],
            centernodes = pipelines,
            outnodes = [fin_node],
            edges = [g.edge(fout_node,pipelines[n],n,0) for n in std.range(0,npipe-1)] + 
            [g.edge(pipelines[n],fin_node,0,n) for n in std.range(0,npipe-1)],
        ),
    }.ret,

    // similar as multifanpipe but jusnt fan-out then pipelines with ending sinks
    multifanout :: function( fout, pipelines,
    fout_nnodes=[1,8,16], fout_multi=[8,2,7],
    name='multifanout', tag_rules=[] ) {
        local fout_nlayers = std.length(fout_multi),
        assert fout_nlayers >= 2 : "fout_nlayers should be >= 2",
        local npipe = std.length(pipelines),
        assert npipe == fout_nnodes[fout_nlayers-1]*fout_multi[fout_nlayers-1] :
            "fout layout error npipe=%d, "%npipe + "fout=%d"%(fout_nnodes[fout_nlayers-1]*fout_multi[fout_nlayers-1]),

        // function to create nodes for one layer
        local fout_layer(ilayer,nnodes,nmulti) = {
            ret : [
                g.pnode({
                    type: fout,
                    name: name+"_fout_%d"%ilayer + "_%d"%inode,
                    data: {
                        multiplicity: nmulti,
                        tag_rules: [],
                    }}, nin=1, nout=nmulti
                ) for inode in std.range(0,nnodes-1)],
        }.ret,
        // nodes for all layers
        local fout_layers = [
            fout_layer(ilayer,
            fout_nnodes[ilayer],
            fout_multi[ilayer])
            for ilayer in std.range(0,fout_nlayers-1)
        ],
        // make edges to make a combo node
        local fout_node = g.intern(
            innodes = fout_layers[0],
            centernodes = if fout_nlayers == 2 then [] else std.flattenArrays([fout_layers[i] for i in std.range(1,fout_nlayers-2)]),
            outnodes = fout_layers[fout_nlayers-1],
            edges = std.flattenArrays(
                [
                    [
                        g.edge(
                            fout_layers[ilayer-1][std.floor(inode/fout_multi[ilayer-1])],
                            fout_layers[ilayer][inode],
                            inode%fout_multi[ilayer-1],
                            0)
                    for inode in std.range(0,fout_nnodes[ilayer]-1)]
                for ilayer in std.range(1,fout_nlayers-1)])
        ),

        // connect comb_fan_out-piples
        ret : g.intern(
            innodes = [fout_node],
            centernodes = pipelines,
            outnodes = [],
            edges = [g.edge(fout_node,pipelines[n],n,0) for n in std.range(0,npipe-1)]
        ),
    }.ret,

}
