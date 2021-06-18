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
}
