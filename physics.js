/* ============================================================
   Bottle-flip physics core (planar rigid body + lumped water)
   ------------------------------------------------------------
   HAND PHASE — arm swing + wrist flick + finger snap:
     φ : hand angle about the wrist W (0 = hanging straight down, CCW +)
     ψ = φ + Δ0 : bottle angle (0 = upright); hand and bottle swing as one
     The wrist W translates at the arm-swing velocity and stops hard at the
     end of the flick; the wrist torque τ rotates hand+bottle about W.
     At release the fingers snap the cap with impulse J = F_snap·Δt:
       Δω = J·r_cap / I_cm   (spin),   Δv = J / M   (slight pull-back)
   FLIGHT / CONTACT — generalized coordinates:
     r = (x, y)  : bottle bottom-centre point (world, m)
     θ           : bottle tilt from upright, CCW positive (rad)
     z_i         : axial position of water lump i from the bottle bottom
   Unit vectors:  e = (-sinθ, cosθ) bottom→cap,  n = de/dθ = (-cosθ, -sinθ)
   Equations of motion (Lagrangian, exact for this model):
     m_b·A                     = F·e − ΣQ_i + θ̇²·m_b·z_b
     M·B  + S·θ̈               = F·n − 2θ̇·Σm_i ż_i
     S·B  + J·θ̈               = τ   − 2θ̇·Σm_i z_i ż_i
     m_i·z̈_i                   = Q_i + m_i z_i θ̇² − m_i·A
   ============================================================ */
var BottleFlip = (function () {
  'use strict';
  var G = 9.81;
  var DEG = Math.PI / 180;

  var BOTTLE = {
    Ltot: 0.225,     // total length incl. cap (m)
    Hsh: 0.185,      // shoulder height (contact corner)
    r: 0.0295,       // body radius (m) → π r² Hliq ≈ 0.5 L
    rCap: 0.0148,
    mEmpty: 0.022,   // empty PET bottle + cap (kg)
    zb: 0.098,       // empty bottle centre of mass height (m)
    Hliq: 0.185,     // liquid column height when full (m)
    zTop: 0.196,     // furthest the water can travel toward the cap
    capacity: 0.5,   // kg of water when full
    zGrip: 0.215     // pinch point (cap centre) from bottom
  };

  var HAND = {
    len: 0.105,        // wrist pivot → pinch point (m)
    mass: 0.50,        // hand mass (kg) for gravity torque
    inertia: 0.006,    // kg·m² effective inertia of the hand (+ part of the forearm) about the wrist
    phi0: -20 * DEG,   // cocked-back start angle
    phiStop: 135 * DEG,// wrist range-of-motion stop
    kStop: 60,         // N·m/rad stop stiffness
    cStop: 0.6,        // N·m·s stop damping
    tStop: 0.05,       // s: the arm swing comes to a hard stop over this time
    tSnap: 0.012,      // s after the arm starts stopping: fingers snap the cap and let go
    snapTime: 0.005    // s: duration of the finger snap (impulse = snap force × snapTime)
  };

  var NUM = {
    dtHand: 1e-4,
    dt: 2e-4,
    frameEvery: 10,      // record a frame every 10 steps (2 ms)
    kGround: 1.5e4,      // N/m
    zetaGround: 0.35,    // damping ratio → restitution ≈ 0.3
    mu: 0.5,             // friction, PET on wood/cloth
    vEps: 0.02,          // friction regularisation (m/s)
    Kacc: 3e4,           // water lump wall stiffness per kg (1/s²)
    zetaWater: 1.0,      // lump/wall damping ratio (inelastic water)
    gamma: 6,            // viscous drag of water along the wall (1/s)
    N: 8,                // water lumps
    tMax: 3.0,           // s after release
    settleTime: 0.15
  };

  function wrapPi(a) {
    a = (a + Math.PI) % (2 * Math.PI);
    if (a < 0) a += 2 * Math.PI;
    return a - Math.PI;
  }

  function defaults() {
    return {
      fill: 0.30,        // 0..1
      torque: 0.6,       // N·m wrist torque during the swing
      snap: 7.5,         // N finger-snap force at the cap (5 ms)
      flickTime: 0.15,   // s of wrist flick / arm swing
      wristHeight: 0.30, // m above the floor
      armSpeed: 1.2,     // m/s arm swing speed (stops hard at the end of the flick)
      armAngle: 70       // deg above horizontal (0..90; a downward swing would drive the bottle into the floor)
    };
  }

  /**
   * Run one throw.
   * opts.fast   – stop as soon as the outcome is certain (for sweeps)
   * opts.frames – record animation frames
   */
  function simulate(p, opts) {
    opts = opts || {};
    var fast = !!opts.fast;
    var recFrames = opts.frames !== false;
    var B = BOTTLE, H = HAND, Q = NUM;
    var i, k;

    /* ---------- masses ---------- */
    var fill = Math.max(0, Math.min(1, p.fill));
    var mw = fill * B.capacity;
    var N = mw > 0.003 ? Q.N : 0;
    var mi = N ? mw / N : 0;
    var hw = fill * B.Hliq;
    var l = N ? hw / N : 0;
    var Ii = mi * (B.r * B.r / 4 + l * l / 12);
    var mb = B.mEmpty, zb = B.zb;
    var Ib = mb * (B.r * B.r / 2 + B.Ltot * B.Ltot / 12);
    var M = mb + mw;
    var z = new Float64Array(N), zd = new Float64Array(N), Qi = new Float64Array(N), zdd = new Float64Array(N);
    for (i = 0; i < N; i++) { z[i] = l / 2 + i * l; zd[i] = 0; }
    var zMin = l / 2, zMax = B.zTop - l / 2;
    var ks = Q.Kacc * mi;
    var csWall = 2 * Q.zetaWater * Math.sqrt(ks * mi);
    var csPair = 2 * Q.zetaWater * Math.sqrt(ks * mi / 2);

    var S0 = mb * zb, J0 = Ib + mb * zb * zb;
    for (k = 0; k < N; k++) { S0 += mi * z[k]; J0 += Ii + mi * z[k] * z[k]; }
    var c0 = S0 / M;                  // CoM height from bottom (water compact)
    var Icm0 = J0 - M * c0 * c0;

    /* ---------- hand phase: arm swing + wrist flick (hand and bottle move as one), finger snap at release ---------- */
    var Lh = H.len, mh = H.mass, Ih = H.inertia;
    var rb = B.zGrip - c0;            // pinch point → bottle CoM (along the axis toward the bottom)
    var aArm = p.armAngle * DEG;
    var vArm = [p.armSpeed * Math.cos(aArm), p.armSpeed * Math.sin(aArm)];
    var W0 = [0, p.wristHeight];
    var T = p.flickTime, ts = H.tStop;
    // wrist trajectory: constant velocity during the flick, then a hard stop over ts
    function wristAt(tt) {
      var sfac, vfac, afac;
      if (tt <= T) { sfac = tt; vfac = 1; afac = 0; }
      else if (tt <= T + ts) { var u = tt - T; sfac = T + u - u * u / (2 * ts); vfac = 1 - u / ts; afac = -1 / ts; }
      else { sfac = T + ts / 2; vfac = 0; afac = 0; }
      return [W0[0] + vArm[0] * sfac, W0[1] + vArm[1] * sfac, vArm[0] * vfac, vArm[1] * vfac, vArm[0] * afac, vArm[1] * afac];
    }
    var kg = Q.kGround, cgH = 2 * Q.zetaGround * Math.sqrt(kg * M);

    var phi = H.phi0, phid = 0, t = 0;
    // the bottle hangs from the fingertips; if it would go through the floor, lean it (bottom corner on the floor)
    var Gy0 = W0[1] - Lh * Math.cos(phi);
    var psi0 = 0, Rc = Math.hypot(B.zGrip, B.r), dc = Math.atan2(B.r, B.zGrip);
    if (Gy0 - 0.001 < B.zGrip) psi0 = dc + Math.acos(Math.max(-1, Math.min(1, (Gy0 - 0.001) / Rc)));
    if (psi0 < phi) psi0 = phi;
    var D0 = psi0 - phi;              // fixed angle between hand and bottle during the swing
    var psi = psi0, psid = 0;
    // inertia of hand + bottle about the wrist
    var qx0 = Lh * Math.sin(phi) + rb * Math.sin(psi), qy0 = -Lh * Math.cos(phi) - rb * Math.cos(psi);
    var rho = Math.hypot(qx0, qy0);
    var Itot = Ih + Icm0 + M * rho * rho;

    var frames = [];
    var handSteps = 0, frameEveryHand = Q.frameEvery * 2;
    var released = false, releaseReason = '';
    var Wc = wristAt(0), hingeF = 0;

    function handFrame() {
      var Gx = Wc[0] + Lh * Math.sin(phi), Gy = Wc[1] - Lh * Math.cos(phi);
      var zz = new Float64Array(N); for (var q = 0; q < N; q++) zz[q] = z[q];
      return { t: t, phase: 0, phi: phi, phid: phid, Wx: Wc[0], Wy: Wc[1],
               x: Gx + B.zGrip * Math.sin(psi), y: Gy - B.zGrip * Math.cos(psi), th: psi, om: psid,
               z: zz, I: Icm0, tilt: wrapPi(psi), contact: false, hinge: hingeF };
    }
    if (recFrames) frames.push(handFrame());

    var dth = Q.dtHand, tRelTarget = T + H.tSnap;
    while (t < tRelTarget - 1e-9) {
      Wc = wristAt(t);
      var sph = Math.sin(phi), cph = Math.cos(phi), sps = Math.sin(psi), cps = Math.cos(psi);
      var ax = -Wc[4], ay = -G - Wc[5];          // effective field in the wrist frame
      var tauW = p.torque;
      if (phi > H.phiStop) tauW += -H.kStop * (phi - H.phiStop) - H.cStop * phid;
      // generalized force of the field on hand CoM and bottle CoM
      var dqx = Lh * cph + rb * cps, dqy = Lh * sph + rb * sps;      // d(CoM)/dφ
      var Qf = (mh * Lh / 2) * (ax * cph + ay * sph) + M * (ax * dqx + ay * dqy);
      // floor contact on the bottle bottom corners
      var Gx = Wc[0] + Lh * sph, Gy = Wc[1] - Lh * cph;
      var Qc = 0, Ffx = 0, Ffy = 0;
      for (var sgn = -1; sgn <= 1; sgn += 2) {
        var sc = sgn * B.r;
        var py = Gy - B.zGrip * cps - sc * sps;
        if (py < 0) {
          var dpx = Lh * cph + B.zGrip * cps + sc * sps, dpy = Lh * sph + B.zGrip * sps - sc * cps;
          var vpx = Wc[2] + phid * dpx, vpy = Wc[3] + phid * dpy;
          var Fn = kg * (-py) - cgH * vpy;
          if (Fn > 0) {
            var sl = vpx / Q.vEps; if (sl > 1) sl = 1; else if (sl < -1) sl = -1;
            var Ft = -Q.mu * Fn * sl;
            Qc += Ft * dpx + Fn * dpy; Ffx += Ft; Ffy += Fn;
          }
        }
      }
      var phidd = (tauW + Qf + Qc) / Itot;
      // finger (hinge) force on the bottle = M·a_c(world) − M·g − F_floor
      var acx = Wc[4] + phidd * dqx - phid * phid * (Lh * sph + rb * sps);
      var acy = Wc[5] + phidd * dqy + phid * phid * (Lh * cph + rb * cps);
      hingeF = Math.hypot(M * acx - Ffx, M * acy + M * G - Ffy);

      phid += phidd * dth; phi += phid * dth; psi = phi + D0; psid = phid;
      t += dth; handSteps++;
      if (recFrames && handSteps % frameEveryHand === 0) frames.push(handFrame());
    }
    releaseReason = 'snap';
    if (recFrames) frames.push(handFrame());

    /* ---------- release: finger snap impulse at the cap ---------- */
    Wc = wristAt(t);
    var tRel = t, phiRel = phi, phidRel = phid, thRel = psi;
    var sph2 = Math.sin(phi), cph2 = Math.cos(phi), sps2 = Math.sin(psi), cps2 = Math.cos(psi);
    var Gx2 = Wc[0] + Lh * sph2, Gy2 = Wc[1] - Lh * cph2;
    var vGx = Wc[2] + phid * Lh * cph2, vGy = Wc[3] + phid * Lh * sph2;
    var om = psid, th = psi;
    var vcx = vGx + om * rb * cps2, vcy = vGy + om * rb * sps2;     // CoM velocity from the swing
    var Jsnap = Math.max(0, p.snap || 0) * H.snapTime;              // impulse (N·s) delivered by the fingers at the cap
    var dOm = Jsnap * rb / Icm0, dV = Jsnap / M;                    // spin up, slight pull-back of the CoM
    om += dOm; vcx -= dV * cps2; vcy -= dV * sps2;
    var swingOmega = psid, snapOmega = dOm;
    // convert CoM state to bottom-point state
    var x = Gx2 + B.zGrip * sps2, y = Gy2 - B.zGrip * cps2;
    var vx = vcx - om * c0 * (-cps2), vy = vcy - om * c0 * (-sps2);  // r = CoM − c0·e, ṙ = v_c − ω·c0·n
    var omRel = om, IRel = Icm0;
    var releaseSpeed = Math.hypot(vcx, vcy);
    var releaseAngle = Math.atan2(vcy, vcx) / DEG;
    var comRelY = Gy2 - rb * cps2;
    var wristAtRelease = wristAt;

    /* ---------- free flight + contact ---------- */
    var dt = Q.dt;
    var cg = cgH;
    var contacts = [[0, B.r], [0, -B.r], [B.Hsh, B.r], [B.Hsh, -B.r], [B.Ltot, B.rCap], [B.Ltot, -B.rCap]];
    var tLand = -1, tiltLand = 0, omLand = 0, ILand = 0, thLand = 0, spreadLand = 0, vLand = 0, vxLand = 0, xLand = 0;
    var apex = comRelY, settledFor = 0, outcome = 'timeout', settled = false;
    var step = 0, tf = 0, everContact = false;
    var maxTiltAfterLand = 0, cEnd = c0;

    while (tf < Q.tMax) {
      var sn = Math.sin(th), cs = Math.cos(th);
      var ex = -sn, ey = cs, nx = -cs, ny = -sn;
      var S = mb * zb, J = Ib + mb * zb * zb, Sd = 0, Szz = 0;
      for (i = 0; i < N; i++) { S += mi * z[i]; J += Ii + mi * z[i] * z[i]; Sd += mi * zd[i]; Szz += mi * z[i] * zd[i]; }

      var Fx = 0, Fy = -M * G, tau = S * G * sn;
      var inContact = false;
      for (k = 0; k < 6; k++) {
        var zk = contacts[k][0], sk = contacts[k][1];
        var pyk = y + zk * ey + sk * ny;
        if (pyk < 0) {
          var vpx2 = vx + om * (zk * nx - sk * ex), vpy2 = vy + om * (zk * ny - sk * ey);
          var Fn2 = kg * (-pyk) - cg * vpy2;
          if (Fn2 > 0) {
            inContact = true;
            var sl2 = vpx2 / Q.vEps; if (sl2 > 1) sl2 = 1; else if (sl2 < -1) sl2 = -1;
            var Ft2 = -Q.mu * Fn2 * sl2;
            Fx += Ft2; Fy += Fn2;
            var rx = zk * ex + sk * nx, ry = zk * ey + sk * ny;
            tau += rx * Fn2 - ry * Ft2;
          }
        }
      }

      var sumQ = 0;
      for (i = 0; i < N; i++) Qi[i] = -mi * G * cs - Q.gamma * mi * zd[i];
      if (N) {
        var pen = zMin - z[0], f;
        if (pen > 0) { f = ks * pen - csWall * zd[0]; if (f < 0) f = 0; Qi[0] += f; }
        pen = z[N - 1] - zMax;
        if (pen > 0) { f = ks * pen + csWall * zd[N - 1]; if (f < 0) f = 0; Qi[N - 1] -= f; }
        for (i = 0; i < N - 1; i++) {
          pen = l - (z[i + 1] - z[i]);
          if (pen > 0) { f = ks * pen - csPair * (zd[i + 1] - zd[i]); if (f < 0) f = 0; Qi[i + 1] += f; Qi[i] -= f; }
        }
        for (i = 0; i < N; i++) sumQ += Qi[i];
      }

      var Fe = Fx * ex + Fy * ey, Fnn = Fx * nx + Fy * ny;
      var A = (Fe - sumQ + om * om * mb * zb) / mb;
      var Fn3 = Fnn - 2 * om * Sd, T2 = tau - 2 * om * Szz;
      var det2 = M * J - S * S;
      var Bacc = (J * Fn3 - S * T2) / det2;
      var thdd = (M * T2 - S * Fn3) / det2;
      for (i = 0; i < N; i++) zdd[i] = Qi[i] / mi + z[i] * om * om - A;

      vx += (A * ex + Bacc * nx) * dt; vy += (A * ey + Bacc * ny) * dt;
      om += thdd * dt;
      for (i = 0; i < N; i++) zd[i] += zdd[i] * dt;
      x += vx * dt; y += vy * dt; th += om * dt;
      for (i = 0; i < N; i++) z[i] += zd[i] * dt;
      tf += dt; t += dt; step++;

      var cy = y + (S / M) * ey; if (cy > apex) apex = cy;
      cEnd = S / M;

      if (inContact && !everContact) {
        everContact = true; tLand = tf; tiltLand = wrapPi(th) / DEG; omLand = om; thLand = th;
        ILand = J - S * S / M; vLand = Math.hypot(vx, vy); vxLand = vx; xLand = x;
        spreadLand = N ? (z[N - 1] - z[0] + l) : 0;
      }
      if (everContact) {
        var tl = Math.abs(wrapPi(th)) / DEG;
        if (tl > maxTiltAfterLand) maxTiltAfterLand = tl;
        var speed = Math.hypot(vx, vy), zq = 0;
        for (i = 0; i < N; i++) if (Math.abs(zd[i]) > zq) zq = Math.abs(zd[i]);
        if (inContact && speed < 0.03 && Math.abs(om) < 0.25 && zq < 0.05) settledFor += dt; else settledFor = 0;
        if (settledFor >= Q.settleTime) { settled = true; outcome = classify(tl); }
        else if (fast && inContact && tl > 70 && tl < 110 && tf - tLand > 0.05) { outcome = 'fallen'; settled = true; }
        else if (tf - tLand > 1.6 && fast) { settled = true; outcome = classify(tl); }
      }

      if (recFrames && step % Q.frameEvery === 0) {
        var zz2 = new Float64Array(N); for (i = 0; i < N; i++) zz2[i] = z[i];
        var Wf = wristAtRelease(t);
        frames.push({ t: t, phase: everContact ? 2 : 1, phi: phiRel, phid: phidRel, Wx: Wf[0], Wy: Wf[1],
                      x: x, y: y, th: th, om: om, z: zz2, I: J - S * S / M, tilt: wrapPi(th), contact: inContact, hinge: 0 });
      }
      if (settled) break;
    }
    function classify(tl) { return tl < 20 ? 'upright' : tl > 160 ? 'inverted' : 'fallen'; }
    if (!settled) { outcome = everContact ? classify(Math.abs(wrapPi(th)) / DEG) : 'fallen'; }

    return {
      params: p, frames: frames, N: N, lumpLen: l, waterHeight: hw, massWater: mw, massTotal: M,
      tRelease: tRel, phiRelease: phiRel, phidRelease: phidRel, releaseReason: releaseReason, thRelease: thRel,
      omegaRelease: omRel, swingOmega: swingOmega, snapOmega: snapOmega, releaseSpeed: releaseSpeed, releaseAngle: releaseAngle, IRelease: IRel, comRelease: c0,
      apex: apex, tLand: tLand, flightTime: tLand > 0 ? tLand : NaN,
      tiltLand: tiltLand, omegaLand: omLand, ILand: ILand, spreadLand: spreadLand, vLand: vLand, vxLand: vxLand,
      xLand: xLand, rotations: tLand > 0 ? (thLand - thRel) / (2 * Math.PI) : NaN,
      maxTiltAfterLand: maxTiltAfterLand, outcome: outcome,
      success: outcome === 'upright', successCap: outcome === 'inverted',
      tEnd: t, tipAngle: Math.atan2(BOTTLE.r, c0) / DEG, tipAngleCap: Math.atan2(BOTTLE.rCap, BOTTLE.Ltot - cEnd) / DEG
    };
  }

  /** Evaluate one cell of the (fill × torque) map. */
  function evalCell(base, fill, snap) {
    var p = {}; for (var k in base) p[k] = base[k];
    p.fill = fill; p.snap = snap;
    var r = simulate(p, { fast: true, frames: false });
    return { fill: fill, snap: snap, tilt: r.tLand > 0 ? Math.abs(r.tiltLand) : 180,
             outcome: r.outcome, success: r.success, successCap: r.successCap, rotations: r.rotations };
  }

  return { simulate: simulate, evalCell: evalCell, defaults: defaults, wrapPi: wrapPi,
           BOTTLE: BOTTLE, HAND: HAND, NUM: NUM, G: G };
})();
