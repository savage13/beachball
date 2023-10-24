"use strict;"

function plot_dc(np1, opts = {}) { //size=200, xy=(0, 0), width=200) {
    /*
    Uses one nodal plane of a double couple to draw a beach ball plot.

    :param ax: axis object of a matplotlib figure
    :param np1: :class:`~NodalPlane`

    Adapted from MATLAB script
    `bb.m <http://www.ceri.memphis.edu/people/olboyd/Software/Software.html>`_
    written by Andy Michael, Chen Ji and Oliver Boyd.
     */
    // check if one or two widths are specified (Circle or Ellipse)
    size = opts.size || 100
    s_1 = np1.strike;
    d_1 = np1.dip;
    r_1 = np1.rake;
    m = 0;
    if (r_1 > 180) {
        r_1 -= 180;
        m = 1;
    }
    if (r_1 < 0) {
        r_1 += 180;
        m = 1;
    }
    //console.log(np1, s_1, d_1, r_1)

    // Get azimuth and dip of second plane
    let aux = aux_plane(s_1, d_1, r_1);
    let s_2 = aux.strike
    let d_2 = aux.dip
    let r_2 = aux.rake;
    d = size / 2;

    if (d_1 >= 90)
        d_1 = 89.9999;
    if (d_2 >= 90)
        d_2 = 89.9999;

    // arange checked for numerical stability, Math.PI is not multiple of 0.1
    phis = arange(0, Math.PI, .01);
    let l1 = []
    let l2 = []
    for(const phi of phis) {
        l1.push( Math.sqrt(
            Math.pow(90 - d_1, 2) / (
                Math.pow(Math.sin(phi), 2) +
                    Math.pow(Math.cos(phi), 2) *
                    Math.pow(90 - d_1, 2) / Math.pow(90, 2))))
        l2.push( Math.sqrt(
            Math.pow(90 - d_2, 2) / (
                Math.pow(Math.sin(phi), 2) + Math.pow(Math.cos(phi), 2) *
                    Math.pow(90 - d_2, 2) / Math.pow(90, 2))))
    }
    //console.log(s_1, d_1, r_1)
    //console.log(s_2, d_2, r_2)
    //console.log(phis, l1, l2)

    const D2R = Math.PI / 180.0;
    let x = []
    let y = []
    let collect = [];
    // plot paths, once for tension areas and once for pressure areas
    for (m_ in [(m + 1) % 2, m]) {
        //console.log('m',m_)
        inc = opts.inc || 10
        tmp = phis.map(p => p + s_1 * D2R)
        ret = pol2cart(tmp, l1);
        x_1 = ret.x
        y_1 = ret.y
        if (m_ == 1) {
            lo = s_1 - 180;
            hi = s_2;
            if (lo > hi)
                inc = -inc;
            //console.log(m_, lo, hi, inc)
            th1 = arange(s_1 - 180, s_2, inc)
            th1.push(s_2)
            //console.log(th1)
            tmp = th1.map(t => t * D2R)
            tmp2 = ones(th1.length).map(v => v * 90)
            ret = pol2cart(tmp, tmp2)
            xs_1 = ret.x
            ys_1 = ret.y
            tmp = phis.map(v => v + s_2 * D2R)
            ret = pol2cart(tmp, l2);
            x_2 = ret.x
            y_2 = ret.y
            th2 = arange(s_2 + 180, s_1, -inc);
        } else {
            hi = s_1 - 180;
            lo = s_2 - 180;
            if (lo > hi)
                inc = -inc;
            //console.log(m_, hi, lo, -inc)
            th1 = arange(hi, lo, -inc)
            th1.push(lo)
            //console.log(th1)
            tmp = th1.map(v => v * D2R)
            tmp2 = ones(th1.length).map(v => v * 90)
            ret = pol2cart(tmp, tmp2)
            xs_1 = ret.x;
            ys_1 = ret.y
            tmp = phis.map(v => v + s_2 * D2R)
            ret = pol2cart(tmp, l2);
            x_2 = ret.x
            y_2 = ret.y
            //x_2 = x_2[::-1];
            //y_2 = y_2[::-1];
            x_2.reverse()
            y_2.reverse()
            th2 = arange(s_2, s_1, inc)
        }
        tmp = th2.map(v => v * D2R)
        tmp2 = ones(th2.length).map(v => v * 90)
        ret = pol2cart(tmp, tmp2)
        xs_2 = ret.x;
        ys_2 = ret.y;
        //console.log(x_1.length, xs_1.length, x_2.length, xs_2.length)
        let xt = [...x_1, ...xs_1, ...x_2, ...xs_2]
        let yt = [...y_1, ...ys_1, ...y_2, ...ys_2]
        //x = Math.concatenate((x_1, xs_1[0], x_2, xs_2[0]));
        //y = Math.concatenate((y_1, ys_1[0], y_2, ys_2[0]));
        xt = xt.map(v => v * d / 90)
        yt = yt.map(v => v * d / 90);
        collect.push({ x: xt, y: yt })
        //x.push(xt)
        //y.push(yt)
        // calculate resolution
        const res = [];
        //for(let value of width)
        //    res.push(value / float(size))

        // construct the patch
        //collect.push(xy2patch(y, x, res, xy));
    }
    pri = { strike: s_1, dip: d_1, rake: r_1 }
    colors = [ opts.color || 'black', opts.color_aux || 'white' ]
    return {collect, colors, aux, pri}//['b', 'w'], collect;
}

function ones(n) {
    let out = []
    for(let i = 0; i < n; i++) {
        out[i] = 1.0
    }
    return out;
}

function zeros(n) {
    if(Array.isArray(n) && n.length == 2) {
        let out = []
        for(let i = 0; i < n[0]; i++) {
            row = []
            for(let j = 0; j < n[1]; j++) {
                row.push(0)
            }
            out.push(row)
        }
        return out
    }
    let out = []
    for(let i = 0; i < n; i++) {
        out[i] = 0.0;
    }
    return out;
}

function arange(start, stop, step) {
    let out = []
    let x0 = start;
    let i = 0;
    if(step == 0.0)
        return out
    while((step > 0 && start < stop) || (step < 0 && start > stop)) {
        out.push(start);
        i += 1;
        start = x0 + i * step;
    }
    return out
}

function pol2cart(th, r) {
    let x = []
    let y = []
    for(let i = 0; i < th.length; i++) {
        x.push( r[i] * Math.cos(th[i]) )
        y.push( r[i] * Math.sin(th[i]) )
    }
    return {x, y};

}

function aux_plane(s1, d1, r1) {
    /*
    Get Strike and dip of second plane.

    Adapted from MATLAB script
    `bb.m <http://www.ceri.memphis.edu/people/olboyd/Software/Software.html>`_
    written by Andy Michael, Chen Ji and Oliver Boyd.
     */
    r2d = 180 / Math.PI;

    z = (s1 + 90.0) / r2d;
    z2 = d1 / r2d;
    z3 = r1 / r2d;
    //console.log(s1,d1,r1,r2d)
    //console.log(z,z2, z3)

    // slick vector in plane 1
    let sl1 = -Math.cos(z3) * Math.cos(z) - Math.sin(z3) * Math.sin(z) * Math.cos(z2);
    let sl2 = Math.cos(z3) * Math.sin(z) - Math.sin(z3) * Math.cos(z) * Math.cos(z2);
    let sl3 = Math.sin(z3) * Math.sin(z2)
    //console.log(sl1, sl2, sl3)
    const ret = strike_dip(sl2, sl1, sl3);
    let strike = ret.strike;
    let dip = ret.dip;
    //console.log(strike, dip)
    n1 = Math.sin(z) * Math.sin(z2)  // normal vector to plane 1
    n2 = Math.cos(z) * Math.sin(z2);
    h1 = -sl2  // strike vector of plane 2
    h2 = sl1;
    // note h3=0 always so we leave it out
    // n3 = Math.cos(z2)

    z = h1 * n1 + h2 * n2;
    z = z / Math.sqrt(h1 * h1 + h2 * h2);
    // we might get above 1.0 only due to floating point
    // precision. Clip for those cases.
    let float64epsilon = 2.2204460492503131e-16;
    if (1.0 < Math.abs(z) < 1.0 + 100 * float64epsilon)
        z = copysign(1.0, z);
    z = Math.acos(z);
    rake = 0;
    if (sl3 > 0)
        rake = z * r2d;
    if (sl3 <= 0)
        rake = -z * r2d;
    return {strike, dip, rake}
}

function copysign(x1, x2) {
    x1 = Math.abs(x1)
    if(x2 < 0.0)
        x1 *= -1
    return x1;
}

function strike_dip(n, e, u) {
    /*
    Finds strike and dip of plane given normal vector having components n, e,
    and u.

    Adapted from MATLAB script
    `bb.m <http://www.ceri.memphis.edu/people/olboyd/Software/Software.html>`_
    written by Andy Michael, Chen Ji and Oliver Boyd.
     */
    r2d = 180 / Math.PI;
    if (u < 0) {
        n = -n;
        e = -e;
        u = -u;
    }
    strike = Math.atan2(e, n) * r2d;
    strike = strike - 90;
    while (strike >= 360)
        strike = strike - 360;
    while (strike < 0)
        strike = strike + 360;
    let x = Math.sqrt(Math.pow(n, 2) + Math.pow(e, 2));
    //console.log(x, u, r2d)
    dip = Math.atan2(x, u) * r2d;

    return {strike, dip}
}

function argmin(v) {
    let imin = 0
    let vmin = v[imin]
    for(let i = 1; i < v.length; i++) {
        if(v[i] < vmin) {
            vmin = v[i]
            imin = i
        }
    }
    return imin
}
function argmax(v) {
    let imax = 0
    let vmax = v[imax]
    for(let i = 1; i < v.length; i++) {
        if(v[i] > vmax) {
            vmax = v[i]
            imax = i
        }
    }
    return imax
}

function mt2plane(mt) {
    /*
    Calculates a nodal plane of a given moment tensor.

    :param mt: :class:`~MomentTensor`
    :return: :class:`~NodalPlane`

    Adapted from MATLAB script
    `bb.m <http://www.ceri.memphis.edu/people/olboyd/Software/Software.html>`_
    written by Andy Michael, Chen Ji and Oliver Boyd.
    */
    let M = [ [mt.xx, mt.xy, mt.xz],
              [mt.xy, mt.yy, mt.yz],
              [mt.xz, mt.yz, mt.zz] ];
    //console.log(JSON.stringify(M))
    //let ret = Math.linalg.eig(mt.mt);
    let ret = jacobi23(M);
    let d = ret.val;
    let v = ret.vec;
    eigsrt(ret.val, ret.vec, 3)
    //console.log(JSON.stringify(d))
    //console.log(JSON.stringify(v))
    //d = [d[1], d[0], d[2]]
    //v = [[ v[1][1], -v[1][0], -v[1][2]],
    //     [ v[2][1], -v[2][0], -v[2][2]],
    //     [-v[0][1],  v[0][0],  v[0][2]]]
    d = ret.val
    v = ret.vec
    imax = argmax(d);
    imin = argmin(d);
    //console.log('d',JSON.stringify(d), imin, imax)
    //console.log('v', JSON.stringify(v[1]),d[1])
    //console.log('v', JSON.stringify(v[0]),d[0])
    //console.log('v', JSON.stringify(v[2]),d[2])
    //console.log('v',JSON.stringify(v))
    //console.log('min,max',imin,imax)
    let ae = [0,0,0]
    let an = [0,0,0]
    //console.log("-")
    //console.log('vmax',JSON.stringify([[v[0][imax],v[1][imax],v[2][imax]],d[imax]]))
    //console.log('vmin',JSON.stringify([[v[0][imin],v[1][imin],v[2][imin]],d[imin]]))
    //console.log("-")
    for(let i = 0; i < 3; i++) {
        ae[i] = (v[i][imax] + v[i][imin]) / Math.sqrt(2.0);
        an[i] = (v[i][imax] - v[i][imin]) / Math.sqrt(2.0);
    }
    //console.log('ae',JSON.stringify(ae))
    //console.log('an',JSON.stringify(an))
    aer = Math.sqrt(Math.pow(ae[0], 2) + Math.pow(ae[1], 2) + Math.pow(ae[2], 2));
    anr = Math.sqrt(Math.pow(an[0], 2) + Math.pow(an[1], 2) + Math.pow(an[2], 2));
    //console.log('aer, anr',aer,anr)
    ae = ae.map(v => v / aer);
    let an1 = 0
    let ae1 = 0
    if (!anr) {
        an = [Math.nan, Math.nan, Math.nan]
    } else
        an = an.map(v => v / anr)
    if (an[2] <= 0.) {
        //console.log(JSON.stringify(an),JSON.stringify(ae))
        an1 = an.map(v => v)
        ae1 = ae.map(v => v)
    } else {
        an1 = [-an[0],-an[1],-an[2]]
        ae1 = [-ae[0],-ae[1],-ae[2]]
    }
    //console.log('an1', an1)
    //console.log('ae1', ae1)
    const ret2 = tdl(an1, ae1);
    //console.log(ret2)
    let ft = ret2.ft;
    let fd = ret2.fd;
    let fl = ret2.fl;
    //console.log(ft, fd, fl)
    return { strike: 360 - ft, dip: fd, rake: 180 - fl};
}


function tdl(an, bn) {
    /*
    Helper function for mt2plane.

    Adapted from MATLAB script
    `bb.m <http://www.ceri.memphis.edu/people/olboyd/Software/Software.html>`_
    written by Andy Michael, Chen Ji and Oliver Boyd.
    */
    let xn = an[0]
    let yn = an[1]
    let zn = an[2]
    let xe = bn[0]
    let ye = bn[1]
    let ze = bn[2]
    let aaa = 1.0 / (1000000)
    let con = 57.2957795
    //console.log(xn,yn,zn, xe,ye,ze)
    if (Math.abs(zn) < aaa) {
        fd = 90.
        axn = Math.abs(xn)
        if (axn > 1.0)
            axn = 1.0
        ft = Math.asin(axn) * con
        st = -xn
        ct = yn
        if (st >= 0. && ct < 0)
            ft = 180. - ft
        if (st < 0. && ct <= 0)
            ft = 180. + ft
        if (st < 0. && ct > 0)
            ft = 360. - ft
        fl = Math.asin(Math.abs(ze)) * con
        sl = -ze
        if (Math.abs(xn) < aaa) {
            cl = xe / yn
        } else {
            cl = -ye / xn
        }
        //console.log(sl, cl, fl)
        if (sl >= 0. && cl < 0)
            fl = 180. - fl
        if (sl < 0. && cl <= 0) {
            fl = fl - 180.
        }
        if (sl < 0. && cl > 0) {
            fl = -fl
        }
    } else {
        if (-zn > 1.0)
            zn = -1.0
        fdh = Math.acos(-zn)
        fd = fdh * con
        sd = Math.sin(fdh)
        // print('sd',sd, fdh)
        if(sd == 0)
            sd = 1e-6

        st = -xn / sd
        ct = yn / sd
        sx = Math.abs(st)
        if(sx > 1.0)
            sx = 1.0
        ft = Math.asin(sx) * con
        if(st >= 0. && ct < 0)
            ft = 180. - ft
        if(st < 0. && ct <= 0)
            ft = 180. + ft
        if(st < 0. && ct > 0)
            ft = 360. - ft
        sl = -ze / sd
        sx = Math.abs(sl)
        if(sx > 1.0)
            sx = 1.0
        fl = Math.asin(sx) * con
        if(st == 0) {
            cl = xe / ct
        } else {
            xxx = yn * zn * ze / sd / sd + ye
            cl = -sd * xxx / xn
            if(ct == 0)
                cl = ye / st
        }
        //console.log(sl, cl, fl)
        if(sl >= 0. && cl < 0)
            fl = 180. - fl
        if(sl < 0. && cl <= 0) {
            fl = fl - 180.
        }
        if(sl < 0. && cl > 0) {
            fl = -fl
        }
    }
    //console.log('ft,fd,fl',ft,fd,fl)
    return {ft, fd, fl}
}

function eye(n) {
    let A = []
    for(let i = 0; i < n; i++) {
        let row = []
        for(let j = 0; j < n; j++) {
            row.push( (i == j) ? 1 : 0 )
        }
        A.push(row)
    }
    return A
}

function diagonal(A) {
    let d = []
    for(let i = 0; i < A.length; i++) {
        d.push(A[i][i])
    }
    return d
}
function finite(A) {
    for(let i = 0; i < A.length; i++) {
        for(let j = 0; j < A[i].length; j++) {
            if(!isFinite(A[i][j]))
               return false
        }
    }
    return true;
}


function clone(input) {
    let result = [];
    input.forEach((row, i) => {
        result.push([]);
        row.forEach(col => {
            result[i].push(col);
        });
    });
    return result;
}

function identity(n) {
    let result = [];
    for (let i = 0; i < n; i++) {
        result.push([]);
        for (let j = 0; j < n; j++) {
            result[i][j] = (j == i ? 1 : 0);
        }
    }
    return result;
}

function transpose(input) {
    let cln = clone(input);
    let n = input[0].length;
    for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
            cln[i][j] = input[j][i];
        }
    }
    return cln;
}

function multiply(A, B) {
    let n = A[0].length;
    let result = identity(n);
    for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
            let sum = 0;
            for (let k = 0; k < n; k++) {
                sum += A[i][k] * B[k][j]
            }
            result[i][j] = sum;
        }
    }
    return result;
}


function isDiagonal(input, epsilon) {
    let n = input[0].length;
    for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
            if (i == j) continue;
            if (Math.abs(input[i][j]) > epsilon)
                return false;
        }
    }
    return true;
}

function clean (input, epsilon) {
    let n = input[0].length;
    let result = clone(input);
    for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
            if (Math.abs(input[i][j]) < epsilon)
                result[i][j] = 0;
            else
                result[i][j] = input[i][j];
        }
    }
    return result;
}


function iterate(S, D, n) {
    // find the indices of the largest off-diagonal element (in magnitude) from D
    let di; let dj;
    for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
            if (i == j) continue;
            if (di === undefined || dj === undefined || Math.abs(D[i][j]) > Math.abs(D[di][dj])) {
                di = i;
                dj = j;
            }
        }
    }
    // find the rotational angle
    let angle;
    if (D[di][di] === D[dj][dj]) {
        if (D[di][dj] > 0)
            angle = Math.PI/4;
        else
            angle = -Math.PI/4;
    } else {
        angle = 0.5 * Math.atan(2 * D[di][dj] / (D[di][di] - D[dj][dj]));
    }

    // compute S1
    let S1 = identity(n);
    S1[di][di] = Math.cos(angle);
    S1[dj][dj] = S1[di][di];
    S1[di][dj] = -Math.sin(angle);
    S1[dj][di] = -S1[di][dj]
    // set new values
    return {
        S: multiply(S, S1),
        D: multiply(multiply(transpose(S1), D), S1)
    };
};

function jacobi(input, iterations, epsilon = 1e-7) {
    let n = input[0].length;
    let D = clone(input);
    let S = identity(n);
    let i = 0;
    for (; i < iterations; i++) {
        let itr = iterate(S, D, n);
        S = itr.S;
        D = itr.D;
        if (isDiagonal(D, epsilon)) {
            D = clean(D, epsilon);
            S = clean(S, epsilon);
            break;
        }
    }
    const it = i
    const val = D.map((v,i) => v[i])

    return { vec: S, val: val, it }
}

function xcircle(center, radius, inc = 10) {
    let x = []
    let y = []
    for(let i = 0; i <= 360; i += inc) {
        x.push(center[0] + radius * Math.cos(i * Math.PI/180.0))
        y.push(center[1] + radius * Math.sin(i * Math.PI/180.0))
    }
    return {x, y}
}

function plot_mt(T, N, P, opts = {} ) {
    // size=200, plot_zerotrace=true, // # noqa
    //x0=0, y0=0, xy=(0, 0), width=200) {

    /*
    Uses a principal axis T, N and P to draw a beach ball plot.

    :param ax: axis object of a matplotlib figure
    :param T: :class:`~PrincipalAxis`
    :param N: :class:`~PrincipalAxis`
    :param P: :class:`~PrincipalAxis`

    Adapted from ps_tensor / utilmeca.c / `Generic Mapping Tools (GMT)`_.

    .. _`Generic Mapping Tools (GMT)`: https://gmt.soest.hawaii.edu
    */
    // check if one or two widths are specified (Circle or Ellipse)
    const EPSILON = 1e-5
    const D2R = Math.PI / 180.0
    plot_zerotrace = (opts.plot_zerotrace !== undefined) ? opts.plot_zerotrace : true;
    xy = [0,0]

    width = [50, 50]
    collect = []
    colors = []
    //res = [ width[0]/size , width[1]/size ]

    b = 1
    big_iso = 0
    j = 1
    j2 = 0
    j3 = 0
    n = 0
    azi = zeros([3, 2])
    x = zeros(400)
    y = zeros(400)
    x2 = zeros(400)
    y2 = zeros(400)
    x3 = zeros(400)
    y3 = zeros(400)
    xp1 = zeros(800)
    yp1 = zeros(800)
    xp2 = zeros(400)
    yp2 = zeros(400)

    a = zeros(3)
    p = zeros(3)
    v = zeros(3)
    a[0] = T.strike
    a[1] = N.strike
    a[2] = P.strike
    p[0] = T.dip
    p[1] = N.dip
    p[2] = P.dip
    v[0] = T.val
    v[1] = N.val
    v[2] = P.val

    // console.log('a',a)
    // console.log('p',p)
    // console.log('v',v)

    vi = (v[0] + v[1] + v[2]) / 3.;
    for(let i = 0; i < 3; i++)
        v[i] = v[i] - vi

    x0 = 0
    y0 = 0
    size = width[0]
    radius_size = size
    // console.log('radius_size',radius_size)

    if(Math.abs(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) < EPSILON) {
        console.log("Pure Explosion / Implosion")
        // pure implosion-explosion
        if (vi > 0.) {
            cir = xcircle(xy, width[0])
            //console.log(xy, cir, width[0])
            collect.push(cir)
            colors.push('black')
        }
        if( vi < 0. ) {
            cir = xcircle(xy, width[0])
            collect.push(cir)
            colors.push('white')
        }
        return { colors, collect }
    }
    if( Math.abs(v[0]) >= Math.abs(v[2]) ) {
        d = 0
        m = 2
    } else {
        d = 2
        m = 0
    }
    if (plot_zerotrace) {
        vi = 0.
    }

    f   = -v[1] / v[d]
    iso = vi / v[d]
    //console.log('iso, f', iso, f)
    //console.log('m,d', m,d)
    // Cliff Frohlich, Seismological Research letters,
    // Vol 7, Number 1, January-February, 1996
    // Unless the isotropic parameter lies in the range
    // between -1 and 1 - f there will be no nodes whatsoever

    if (iso < -1) {
        //console.log('iso < -1')
        cir = xcircle(xy, width[0])
        collect.push(cir)
        colors.append('white')
        return colors, collect
    } else if (iso > 1 - f ) {
        //console.log('iso < 1 - f')
        cir = xcircle(xy, width[0])
        collect.push(cir)
        colors.push('blacl')
        return {colors, collect }
    }

    spd = Math.sin(p[d] * D2R)
    cpd = Math.cos(p[d] * D2R)
    spb = Math.sin(p[b] * D2R)
    cpb = Math.cos(p[b] * D2R)
    spm = Math.sin(p[m] * D2R)
    cpm = Math.cos(p[m] * D2R)
    sad = Math.sin(a[d] * D2R)
    cad = Math.cos(a[d] * D2R)
    sab = Math.sin(a[b] * D2R)
    cab = Math.cos(a[b] * D2R)
    sam = Math.sin(a[m] * D2R)
    cam = Math.cos(a[m] * D2R)

    for(let i = 0; i < 360; i++) {
        fir = i * D2R
        s2alphan = (2. + 2. * iso) / (3. + (1. - 2. * f) * Math.cos(2. * fir))
        if(s2alphan > 1.) {
            big_iso += 1
        } else {
            alphan = Math.asin(Math.sqrt(s2alphan))
            sfi = Math.sin(fir)
            cfi = Math.cos(fir)
            san = Math.sin(alphan)
            can = Math.cos(alphan)

            xz = can * spd + san * sfi * spb + san * cfi * spm
            xn = can * cpd * cad + san * sfi * cpb * cab + san * cfi * cpm * cam
            xe = can * cpd * sad + san * sfi * cpb * sab + san * cfi * cpm * sam

            takeoff = 0.0
            if (Math.abs(xn) < EPSILON && Math.abs(xe) < EPSILON) {
                takeoff = 0.
                az = 0.
            } else {
                az = Math.atan2(xe, xn)
                if(az < 0.) {
                    az += Math.PI * 2.0
                    takeoff = Math.acos(xz / (Math.sqrt(xz * xz + xn * xn + xe * xe)))
                }
            }
            if(takeoff > Math.PI / 2.) {
                takeoff = Math.PI - takeoff
                az += Math.PI
                if (az > Math.PI * 2.)
                    az -= Math.PI * 2.
            }
            r = Math.sqrt(2) * Math.sin(takeoff / 2.)
            si = Math.sin(az)
            co = Math.cos(az)
            if (i == 0) {
                azi[i][0] = az
                x[i] = x0 + radius_size * r * si
                y[i] = y0 + radius_size * r * co
                azp = az
            } else {
                //console.log('az-azp',i, Math.abs(az-azp)-Math.PI)
                if (Math.abs(Math.abs(az - azp) - Math.PI) < D2R * 10.) {
                    azi[n][1] = azp
                    n += 1
                    //console.log('n',n)
                    azi[n][0] = az
                }
                if (Math.abs(Math.abs(az - azp) - Math.PI * 2.) < D2R * 2.) {
                    if (azp < az) {
                        azi[n][0] += Math.PI * 2.
                    } else {
                        azi[n][0] -= Math.PI * 2.
                    }
                }
                if (n == 0) {
                    x[j] = x0 + radius_size * r * si
                    y[j] = y0 + radius_size * r * co
                    j += 1
                } else if(n == 1) {
                    x2[j2] = x0 + radius_size * r * si
                    y2[j2] = y0 + radius_size * r * co
                    j2 += 1
                } else if(n == 2) {
                    x3[j3] = x0 + radius_size * r * si
                    y3[j3] = y0 + radius_size * r * co
                    j3 += 1
                }
                azp = az
            }
        }
    }
    azi[n][1] = az
    //console.log(n)
    if (v[1] < 0.) {
        rgb1 = 'black'
        rgb2 = 'white'
    } else {
        rgb1 = 'white'
        rgb2 = 'black'
    }

    cir = xcircle(xy, width[0])
    collect.push(cir)
    colors.push(rgb2)
    if(n == 0) {
        collect.push( {x, y} )
        colors.push(rgb1)
        return { colors, collect }
    } else if(n == 1) {
        for(let i = 0; i < j; i++) {
            xp1[i] = x[i]
            yp1[i] = y[i]
        }
        i = j
        if(azi[0][0] - azi[0][1] > Math.PI) {
            azi[0][0] -= Math.PI * 2.
        } else if(azi[0][1] - azi[0][0] > Math.PI) {
            azi[0][0] += Math.PI * 2.
        }
        if(azi[0][0] < azi[0][1]) {
            az = azi[0][1] - D2R
            while(az > azi[0][0]) {
                si = Math.sin(az)
                co = Math.cos(az)
                xp1[i] = x0 + radius_size * si
                yp1[i] = y0 + radius_size * co
                i += 1
                az -= D2R
            }
        } else {
            az = azi[0][1] + D2R
            while(az < azi[0][0]) {
                si = Math.sin(az)
                co = Math.cos(az)
                xp1[i] = x0 + radius_size * si
                yp1[i] = y0 + radius_size * co
                i += 1
                az += D2R
            }
        }
        collect.push( {x: xp1.slice(0,i), y: yp1.slice(0,i) })
        colors.push(rgb1)
        for(let i = 0; i < j2; i++) {
            xp2[i] = x2[i]
            yp2[i] = y2[i]
        }
        i = j2-1
        //console.log('j2',i,j2)
        if(azi[1][0] - azi[1][1] > Math.PI) {
            azi[1][0] -= Math.PI * 2.
        } else if(azi[1][1] - azi[1][0] > Math.PI) {
            azi[1][0] += Math.PI * 2.
        }
        if(azi[1][0] < azi[1][1]) {
            az = azi[1][1] - D2R
            while(az > azi[1][0]) {
                si = Math.sin(az)
                co = Math.cos(az)
                xp2[i] = x0 + radius_size * si
                i += 1
                yp2[i] = y0 + radius_size * co
                az -= D2R
            }
        } else {
            az = azi[1][1] + D2R
            while(az < azi[1][0]) {
                si = Math.sin(az)
                co = Math.cos(az)
                xp2[i] = x0 + radius_size * si
                i += 1
                yp2[i] = y0 + radius_size * co
                az += D2R
            }
        }
        collect.push({ x: xp2.slice(0,i), y: yp2.slice(0,i) })
        colors.push(rgb1)
        return { colors, collect }
    } else if(n == 2) {
        for(let i  = 0; i < j3; i++) {
            xp1[i] = x3[i]
            yp1[i] = y3[i]
        }
        i = j3
        for(let ii = 0; ii < j; ii++) {
            xp1[i] = x[ii]
            i += 1
            yp1[i] = y[ii]
        }
        if(big_iso) {
            ii = j2 - 1
            while(ii >= 0) {
                xp1[i] = x2[ii]
                i += 1
                yp1[i] = y2[ii]
                ii -= 1
            }
            collect.push({x: xp1.slice(0,i), y: yp1.slice(0,i) })
            colors.push(rgb1)
            return { colors, collect }
        }
        if(azi[2][0] - azi[0][1] > Math.PI) {
            azi[2][0] -= Math.PI * 2.
        } else if(azi[0][1] - azi[2][0] > Math.PI) {
            azi[2][0] += Math.PI * 2.
        }
        if(azi[2][0] < azi[0][1]) {
            az = azi[0][1] - D2R
            while(az > azi[2][0]) {
                si = Math.sin(az)
                co = Math.cos(az)
                xp1[i] = x0 + radius_size * si
                i += 1
                yp1[i] = y0 + radius_size * co
                az -= D2R
            }
        } else {
            az = azi[0][1] + D2R
            while(az < azi[2][0]) {
                si = Math.sin(az)
                co = Math.cos(az)
                xp1[i] = x0 + radius_size * si
                i += 1
                yp1[i] = y0 + radius_size * co
                az += D2R
            }
        }
        collect.push({x: xp1.slice(0,i), y: yp1.slice(0,i) })
        colors.push(rgb1)

        for(let  i = 0; i < j2; i++) {
            xp2[i] = x2[i]
            yp2[i] = y2[i]
        }
        if(azi[1][0] - azi[1][1] > Math.PI) {
            azi[1][0] -= Math.PI * 2.
        } else if(azi[1][1] - azi[1][0] > Math.PI) {
            azi[1][0] += Math.PI * 2.
        }
        if(azi[1][0] < azi[1][1]) {
            az = azi[1][1] - D2R
            while(az > azi[1][0]) {
                si = Math.sin(az)
                co = Math.cos(az)
                xp2[i] = x0 + radius_size * si
                i += 1
                yp2[i] = y0 + radius_size * co
                az -= D2R
            }
        } else {
            az = azi[1][1] + D2R
            while(az < azi[1][0]) {
                si = Math.sin(az)
                co = Math.cos(az)
                xp2[i] = x0 + radius_size * si
                i += 1
                yp2[i] = y0 + radius_size * co
                az += D2R
            }
        }
        collect.push({ x: xp2.slice(0,i), y: yp2.slice(0,i) })
        colors.push(rgb1)
        return { colors, collect }
    }
}

function eigsrt(d, v, n, order = 1){
    let k = 1
    for (let i=1;i<n;i++) {
        let p = d[i-1];
        k = i
        for (let j=i+1;j<=n;j++) {
            if(order > 0) {
                if (d[j-1] <= p) {
                    p = d[j-1];
                    k = j
                }
            }
            if(order < 0) {
                if (d[j-1] >= p) {
                    p = d[j-1];
                    k = j
                }
            }
        }
        if (k != i) {
            d[k-1] = d[i-1];
            d[i-1] = p;
            for (let j=1;j<=n;j++) {
                p = v[j-1][i-1];
                v[j-1][i-1] = v[j-1][k-1];
                v[j-1][k-1] = p;
            }
        }
    }
}


function jacobi23(a) {
    let n = 3
    let b = [0,0,0];
    let z = [0,0,0];
    let d = [0,0,0];
    let v = [[0,0,0],[0,0,0],[0,0,0]]
    for (let ip=1;ip<=n;ip++) {
        for (let iq=1;iq<=n;iq++) {
            v[ip-1][iq-1] = 0.0;
        }
        v[ip-1][ip-1] = 1.0;
    }
    for (let ip=1;ip<=n;ip++) {
        b[ip-1]=d[ip-1] = a[ip-1][ip-1];
        z[ip-1]=0.0;
    }
    let nrot=0;
    for (let i=1;i<=50;i++) {
        //console.log('iter',i, JSON.stringify(v))
        let sm=0.0;
        for (let ip=1;ip<=n-1;ip++) {
            for (let iq=ip+1;iq<=n;iq++)
                sm += Math.abs(a[ip-1][iq-1]);
        }
        if (sm == 0.0) {
            //free_vector(z,1,n);
            //console.log('return',JSON.stringify(v),d)
            //free_vector(b,1,n);
            return { vec: v, val: d }
        }
        let tresh = 0.0
        if (i < 4) {
            tresh=0.2*sm/(n*n);
        }
        // console.log('iter tresh', tresh)
        for (let ip=1;ip<=n-1;ip++) {
            for (let iq=ip+1;iq<=n;iq++) {
                let g=100.0*Math.abs(a[ip-1][iq-1]);
                //console.log('iter ip,iq,g',ip,iq,g)
                if (i > 4 && (Math.abs(d[ip-1])+g) == Math.abs(d[ip-1])
                    && (Math.abs(d[iq-1])+g) == Math.abs(d[iq-1])) {
                    //console.log("iter set = 0", ip, iq)
                    a[ip-1][iq-1]=0.0;
                } else if (Math.abs(a[ip-1][iq-1]) > tresh) {
                    //console.log("iter a > tresh", ip, iq, a[ip-1][iq-1], '>', tresh)
                    let h=d[iq-1]-d[ip-1];
                    let t = 0
                    if ((Math.abs(h)+g) == Math.abs(h)) {
                        t=(a[ip-1][iq-1])/h;
                    } else {
                        let theta=0.5*h/(a[ip-1][iq-1]);
                        t=1.0/(Math.abs(theta)+Math.sqrt(1.0+theta*theta));
                        if (theta < 0.0) {
                            t = -t;
                        }
                    }
                    //console.log("iter t",t);
                    let c=1.0/Math.sqrt(1+t*t);
                    let s=t*c;
                    let tau=s/(1.0+c);
                    h=t*a[ip-1][iq-1];
                    //console.log("iter h",h,s,tau)
                    z[ip-1] -= h;
                    z[iq-1] += h;
                    d[ip-1] -= h;
                    d[iq-1] += h;
                    a[ip-1][iq-1]=0.0;
                    for (let j=1;j<=ip-1;j++) {
                        ROTATE(a,j,ip,j,iq,s,tau)
                    }
                    for (let j=ip+1;j<=iq-1;j++) {
                        ROTATE(a,ip,j,j,iq,s,tau)
                    }
                    for (let j=iq+1;j<=n;j++) {
                        ROTATE(a,ip,j,iq,j,s,tau)
                    }
                    for (let j=1;j<=n;j++) {
                        //console.log("iter ", i, j, ip, iq)
                        ROTATE(v,j,ip,j,iq,s,tau)
                    }
                    //console.log('iter', i, JSON.stringify(v))
                    nrot += 1
                }
            }
        }
        for (let ip=1;ip<=n;ip++) {
            b[ip-1] += z[ip-1];
            d[ip-1] = b[ip-1];
            z[ip-1] = 0.0;
        }
    }
    nrerror("Too many iterations in routine jacobi");
    return {}
}

function ROTATE(a,i,j,k,l, s, tau) {
    let g = a[i-1][j-1]
    let h = a[k-1][l-1]
    a[i-1][j-1] = g - s * (h + g * tau)
    a[k-1][l-1] = h + s * (g - h * tau)
}

function mt2axes(mt) {
    /*
    Calculates the principal axes of a given moment tensor.

    :param mt: :class:`~MomentTensor`
    :return: tuple of :class:`~PrincipalAxis` T, N and P

    Adapted from ps_tensor / utilmeca.c /
    `Generic Mapping Tools (GMT) <https://gmt.soest.hawaii.edu>`_.
    */
    const np = 3
    const R2D = 180.0 / Math.PI;
    let M = [ [mt.xx, mt.xy, mt.xz],
              [mt.xy, mt.yy, mt.yz],
              [mt.xz, mt.yz, mt.zz] ];
    //(d, v) = np.linalg.eigh(M)
    //console.log('M',JSON.stringify(M))
    //console.log('M',M)
    let ret = jacobi23(M)
    eigsrt(ret.val, ret.vec, 3, -1)
    //console.log(ret.val)
    //console.log(transpose(ret.vec))
    //console.log('ret',ret)
    //ret.vec = transpose(ret.vec)
    let d = ret.val
    let v = ret.vec

    let pl = [0,0,0]
    let az = [0,0,0]
    for (let j = 0; j < np; j++) {
        pl[j] = Math.asin(-v[0][j]);
        az[j] = Math.atan2(v[2][j], -v[1][j]);
        if (pl[j] <= 0.) {
            pl[j] = -pl[j];
            az[j] += Math.PI;
        }
        if (az[j] < 0)
            az[j] += 2.0 * Math.PI
        else if (az[j] > 2.0 * Math.PI)
            az[j] -= 2.0 * Math.PI
        pl[j] *= R2D;
        az[j] *= R2D;
    }
    let t = { val: d[0], strike: az[0], dip: pl[0] }
    let n = { val: d[1], strike: az[1], dip: pl[1] }
    let p = { val: d[2], strike: az[2], dip: pl[2] }
    //console.log('t',t)
    //console.log('n',n)
    //console.log('p',p)
    return {t, n, p}
}

function ax2dc(t, p) {
    const EPSIL = 1e-5;
    const R2D = 180.0 / Math.PI
    const D2R = Math.PI / 180.0
    //console.log('in',t.strike, t.dip, p.strike, p.dip)

    let pp = p.strike * D2R
    let dp = p.dip * D2R
    let pt = t.strike * D2R
    let dt = t.dip * D2R
    //console.log('in',pp,dp,pt,dt)
    let sdp = Math.sin(dp)
    let cdp = Math.cos(dp)
    let sdt = Math.sin(dt)
    let cdt = Math.cos(dt)
    let spt = Math.sin(pt)
    let cpt = Math.cos(pt)
    let spp = Math.sin(pp)
    let cpp = Math.cos(pp)

    cpt *= cdt
    spt *= cdt
    cpp *= cdp
    spp *= cdp
    //console.log("sd", sdt, sdp, spt, spp, cpt, cpp)
    let amz = sdt + sdp
    let amx = spt + spp
    let amy = cpt + cpp
    //console.log("am", amx, amy, amz)
    let d1 = Math.atan2(Math.sqrt(amx*amx + amy*amy), amz)
    let p1 = 0.0
    if(Math.abs(amy) < EPSIL && Math.abs(amx) < EPSIL) {
        p1 = 0.0
    } else {
        p1 = Math.atan2(amy, -amx)
    }
    //console.log('p1,d1',p1,d1)
    if (d1 > Math.PI/2 && (Math.abs(d1 - Math.PI/2) > EPSIL)) {
        //console.log("FIX UP p1,d1")
        d1 = Math.PI - d1
        p1 += Math.PI
        if (p1 >= 2.0 * Math.PI)
            p1 -= 2.0 * Math.PI
    }
    if (p1 < 0.)
        p1 += 2.0 * Math.PI
    //console.log('p1,d1',p1,d1)
    amz = sdt - sdp
    amx = spt - spp
    amy = cpt - cpp
    //console.log("am", amx, amy, amz)
    let d2 = Math.atan2(Math.sqrt(amx*amx + amy*amy), amz)
    let p2 = Math.atan2(amy, -amx)
    //console.log("p2,d2", p2, d2)
    if (d2 > (Math.PI/2)) {// && Math.abs(d2 - Math.PI/2) > EPSIL) {
        d2 = Math.PI - d2
        p2 += Math.PI
        //console.log("FIX UP p2,d2", p2, 2*Math.PI, p2 - (2.0 * Math.PI))
        if (p2 >= (2.0 * Math.PI))
            p2 -= 2.0 * Math.PI
    }
    if (p2 < 0.)
        p2 += (2.0 * Math.PI)
    //console.log("p2,d2", p2, d2)
    let np1 = { strike: p1 * R2D, dip: d1 * R2D }
    let np2 = { strike: p2 * R2D, dip: d2 * R2D }

    let im = 1
    if (dp > dt && Math.abs(dp-dt) > EPSIL)
        im = -1
    np1.rake = computed_rake2(np2.strike, np2.dip, np1.strike, np1.dip, im)
    np2.rake = computed_rake2(np1.strike, np1.dip, np2.strike, np2.dip, im)
    //console.log("np1",np1.strike.toFixed(2), np1.dip.toFixed(2), np1.rake.toFixed(2), p1.toFixed(2), im,dp * R2D)
    //console.log("np2",np2.strike.toFixed(2), np2.dip.toFixed(2), np2.rake.toFixed(2), p2.toFixed(2), im,dt * R2D)
    return { NP1: np1, NP2: np2 }
}

function computed_rake2(str1, dip1, str2, dip2, fault) {

/*
   Compute rake in the second nodal plane when strike and dip
   for first and second nodal plane are given with a double
   characterizing the fault :
              +1. inverse fault
              -1. normal fault.
   Angles are in degrees.
*/

    /* Genevieve Patau */
    const EPSIL = 1e-5;
    const D2R = Math.PI / 180.0
    let ss = Math.sin((str1 - str2) * D2R)
    let cs = Math.cos((str1 - str2) * D2R)

    let sd = Math.sin(dip1 * D2R);
    let cd = Math.cos(dip2 * D2R);
    let sinrake2;
    if (Math.abs(dip2 - 90.) < EPSIL) {
        //console.log("cr2: dip2 ~ 90")
        //sinrake2 = fault * cd;
        //sinrake2 = -fault * sd * cs / Math.cos((90.0 - EPSIL)*D2R)
        cd = Math.cos((90.0 - EPSIL) * D2R)
    }
    if(Math.abs(dip1) < EPSIL) {
        sd = Math.sin(EPSIL * D2R)
    }
    sinrake2 = -fault * sd * cs / cd;
    //console.log("cr2:", sinrake2, - fault * sd * ss, fault);
    let rake2 = Math.atan2(sinrake2, -fault * sd * ss) * 180.0 / Math.PI;
    //console.log('cr2',sinrake2.toFixed(2), (-fault * sd * ss).toFixed(2), '::',
    //            ss.toFixed(2), cs.toFixed(2), sd.toFixed(2), cd.toFixed(2), '::', rake2, fault)
    return(rake2);
}


function plotmt(M, options) {
    const size = options.size || 100
    const tcolor = options.tcolor || 'black'
    const pcolor = options.pcolor || 'white'
    const EPSILON = 1e-5
    const ax = mt2axes(M)
    //console.log(ax)
    if(Math.abs(ax.n.val) < EPSILON && Math.abs(ax.t.val + ax.p.val) < EPSILON ) {
        const np = ax2dc(ax.t, ax.p);
        //console.log(np.NP1.strike.toFixed(2), np.NP1.dip.toFixed(2), np.NP1.rake.toFixed(2))
        //console.log(np.NP2.strike.toFixed(2), np.NP2.dip.toFixed(2), np.NP2.rake.toFixed(2))
        //console.log("mech")
        return ps_mechanism(0,0,np,size, tcolor,pcolor,1)
    } else {
        //console.log("tensor", ax.n.val, ax.t.val + ax.p.val)
        //return plot_mt(ax.t, ax.n, ax.p)
        return ps_tensor(0,0, size, ax.t, ax.n, ax.p, tcolor, pcolor, true, false );
    }
}

function null_axis_dip(str1, dip1, str2, dip2) {
/*
   compute null axis dip when strike and dip are given
   for each nodal plane.
   Angles are in degrees.
*/

/* Genevieve Patau */
    const D2R = Math.PI / 180.0
    const R2D = 180.0 / Math.PI
    let den = Math.asin(Math.sin(dip1 * D2R) * Math.sin(dip2  * D2R) * Math.sin((str1 - str2)  * D2R)) * R2D
    if (den < 0.0)
        den = -den;
    return den
}
function null_axis_strike(str1, dip1, str2, dip2) {

/*
   Compute null axis strike when strike and dip are given
   for each nodal plane.
   Angles are in degrees.
*/

    /* Genevieve Patau */
    const D2R = Math.PI / 180.0
    const R2D = 180.0 / Math.PI
    //double phn, cosphn, sinphn;
    //double sd1, cd1, sd2, cd2, ss1, cs1, ss2, cs2;

    sd1 = Math.sin(dip1 * D2R)
    cd1 = Math.cos(dip1 * D2R)
    sd2 = Math.sin(dip2 * D2R)
    cd2 = Math.cos(dip2 * D2R)

    ss1 = Math.sin(str1 * D2R)
    cs1 = Math.cos(str1 * D2R)
    ss2 = Math.sin(str2 * D2R)
    cs2 = Math.cos(str2 * D2R)

    cosphn = sd1 * cs1 * cd2 - sd2 * cs2 * cd1;
    sinphn = sd1 * ss1 * cd2 - sd2 * ss2 * cd1;
    if (Math.sin((str1 - str2) * D2R) < 0.) {
        cosphn = -cosphn;
        sinphn = -sinphn;
    }
    phn = Math.atan2(sinphn, cosphn) * R2D
    if (phn < 0.0)
        phn += 360.;
    return phn
}

function proj_radius(str1, dip1, str) {

/*
   Compute the vector radius for a given strike,
   equal area projection, inferior sphere.
   Strike and dip of the plane are given.
*/

/* Genevieve Patau */

    //double dip, r;
    const D2R = Math.PI / 180.0
    const EPSIL = 1e-5

    if (Math.abs(dip1 - 90.) < EPSIL) {
      /*
        printf("\nVertical plane : strike is constant.");
        printf("\nFor ps_mechanism r == 1 for str = str1");
        printf("\n            else r == 0. is used.");
      */
      r = (Math.abs(str - str1) < EPSIL || Math.abs(str - str1 - 180) < EPSIL) ? 1. : 0.;
    } else {
        dip = Math.atan(Math.tan(dip1 * D2R) * Math.sin((str - str1) * D2R));
        r = Math.sqrt(2.) * Math.sin(Math.PI/4.0 - dip/2.);
    }
    return r;
}

function zero_360(str) {
    if(str >= 360)
        str -= 360
    if(str < 0.0)
        str += 360.0
    return str
}

function computed_strike1(np1) {
    const EPSIL = 1e-5;
    const D2R = Math.PI / 180.0
    const R2D = 180.0 / Math.PI
    let cd1 = Math.cos(np1.dip * D2R)
    let am = (Math.abs(np1.rake) < EPSIL) ? 1.0 : np1.rake / Math.abs(np1.rake);

    let sr = Math.sin(np1.rake * D2R)
    let cr = Math.cos(np1.rake * D2R)
    let ss = Math.sin(np1.strike * D2R)
    let cs = Math.cos(np1.strike * D2R)

    if(cd1 < EPSIL && Math.abs(cr) < EPSIL) {
        return np1.strike + 180.0;
    }
    let temp = cr * cs
    temp += sr * ss * cd1
    sp2 = -am * temp
    temp = ss * cr
    temp -= sr * cs * cd1
    cp2 = am * temp
    str2 = Math.atan2(sp2, cp2) * R2D;
    return zero_360(str2)
}

function computed_dip1(np1) {
    const EPSIL = 1e-5
    const R2D = 180.0 / Math.PI
    const D2R = Math.PI / 180.0
    let am = (Math.abs(np1.rake) < EPSIL) ? 1.0 : np1.rake / Math.abs(np1.rake);
    return Math.acos(am * Math.sin(np1.rake * D2R) * Math.sin(np1.dip * D2R) ) * R2D
}

function computed_rake1(np1, strike2, dip2) {
    const EPSIL = 1e-5;
    const D2R = Math.PI / 180.0
    const R2D = 180.0 / Math.PI
    let am = (Math.abs(np1.rake) < EPSIL) ? 1.0 : np1.rake / Math.abs(np1.rake);
    let sd = Math.sin(np1.dip * D2R)
    let cd = Math.cos(np1.dip * D2R)
    let ss = Math.sin((np1.strike - strike2) * D2R)
    let cs = Math.cos((np1.strike - strike2) * D2R)

    let sinrake2;
    if(Math.abs(dip2 - 90) < EPSIL) {
        sinrake2 = am * cd
    } else {
        sinrake2 = -am * sd * cs / cd
    }
    return Math.atan2(sinrake2, -am * sd * ss)
}

function define_second_plane(np1) {
    let strike = computed_strike1(np1)
    let dip = computed_dip1(np1)
    let rake = computed_rake1(np1, strike, dip)
    return { strike, dip, rake }
}

function ps_mechanism(x0, y0, meca, size, rgb, ergb, outline)  {

    const D2R = Math.PI / 180.0
    const EPSIL = 1e-5
    //double x[1000], y[1000];
    //double pos_NP1_NP2 = sind(meca.NP1.str - meca.NP2.str);
    //double fault = (GMT_IS_ZERO (meca.NP1.rake) ? meca.NP2.rake / Math.abs(meca.NP2.rake) : meca.NP1.rake / Math.abs(meca.NP1.rake));
    //double radius_size;
    //double str, radius, increment;
    //double si, co;

    let collect = []
    let color = []

    pos_NP1_NP2 = Math.sin((meca.NP1.strike-meca.NP2.strike)*D2R)

    let fault = 1
    if(Math.abs(meca.NP1.rake) < EPSIL)
        fault = meca.NP2.rake / Math.abs(meca.NP2.rake)
    else
        fault = meca.NP1.rake / Math.abs(meca.NP1.rake)

    //GMT_LONG lineout = 1, i;
    //GMT_LONG npoints;

    //struct AXIS N_axis;
    let N_axis = { strike: 0, dip: 0, rake: 0};
    /* compute null axis strike and dip */
    N_axis.dip = null_axis_dip(meca.NP1.strike, meca.NP1.dip, meca.NP2.strike, meca.NP2.dip);
    if (Math.abs(90. - N_axis.dip) < EPSIL)
        N_axis.strike = meca.NP1.strike;
    else
        N_axis.strike = null_axis_strike(meca.NP1.strike, meca.NP1.dip, meca.NP2.strike, meca.NP2.dip);

    /* compute radius size of the bubble */
    radius_size = size * 0.5;

    /* outline the bubble */
    //ps_plot(x0 + radius_size, y0, 3);
    /*  argument is DIAMETER!!*/
    //ps_circle(x0, y0, radius_size*2., ergb, lineout);
    let x = []
    let y = []
    collect.push( xcircle([x0,y0], radius_size) )
    color.push(ergb)


    if (Math.abs(pos_NP1_NP2) < EPSIL) {
        //console.log("pure compression/extensional fault")
        /* pure normal or inverse fault (null axis strike is determined
           with + or - 180 degrees. */
        /* first nodal plane part */
        i = -1;
        increment = 1.;
        str = meca.NP1.strike;
        while (str <= meca.NP1.strike + 180. + EPSIL) {
            i++;
            radius = proj_radius(meca.NP1.strike, meca.NP1.dip, str) * radius_size;
            //sincosd (str, &si, &co);
            si = Math.sin(str * D2R)
            co = Math.cos(str * D2R)
            x[i] = x0 + radius * si;
            y[i] = y0 + radius * co;
            str += increment;
        }
        if (fault < 0.) {
            /* normal fault, close first compressing part */
            str = meca.NP1.strike + 180.;
            while (str >= meca.NP1.strike - EPSIL) {
                i++;
                //sincosd (str, &si, &co);
                si = Math.sin(str * D2R)
                co = Math.cos(str * D2R)
                x[i] = x0 + si * radius_size;
                y[i] = y0 + co * radius_size;
                str -= increment;
            }
            npoints = i + 1;
            //ps_polygon (x, y, npoints, rgb, outline);
            collect.push({x: x.slice(0,npoints), y: y.slice(0,npoints) })
            color.push(rgb)
            i = -1;
        }
        /* second nodal plane part */
        str = meca.NP2.strike;
        while (str <= meca.NP2.strike + 180. + EPSIL) {
            i++;
            radius = proj_radius(meca.NP2.strike, meca.NP2.dip, str) * radius_size;
            //sincosd (str, &si, &co);
            si = Math.sin(str * D2R)
            co = Math.cos(str * D2R)
            x[i] = x0 + radius * si;
            y[i] = y0 + radius * co;
            str += increment;
        }
        if (fault > 0.) {
            /* inverse fault, close compressing part */
            npoints = i+1;
            //ps_polygon(x, y, npoints, rgb, outline);
            collect.push({x: x.slice(0,npoints), y: y.slice(0,npoints) })
            color.push(rgb)
        }
        else {
            /* normal fault, close second compressing part */
            str = meca.NP2.strike + 180.;
            while (str >= meca.NP2.strike - EPSIL) {
                i++;
                //sincosd (str, &si, &co);
                si = Math.sin(str * D2R)
                co = Math.cos(str * D2R)
                x[i] = x0 + si * radius_size;
                y[i] = y0 + co * radius_size;
                str -= increment;
            }
            npoints = i + 1;
            collect.push({x: x.slice(0,npoints), y: y.slice(0,npoints) })
            color.push(rgb)
            //ps_polygon(x, y, npoints, rgb, outline);
        }
    }
    /* pure strike-slip */
    else if ((90. - meca.NP1.dip) < EPSIL && (90. - meca.NP2.dip) < EPSIL) {
        //console.log("pure strike slip fault")
        increment = Math.abs(meca.NP1.rake) < EPSIL ? 1. : -1.;
        /* first compressing part */
        i = 0;
        str = meca.NP1.strike;
        while (increment > 0. ? str <= meca.NP1.strike + 90. : str >= meca.NP1.strike - 90.) {
            //sincosd (str, &si, &co);
            si = Math.sin(str * D2R)
            co = Math.cos(str * D2R)
            x[i] = x0 + si * radius_size;
            y[i] = y0 + co * radius_size;
            str += increment;
            i++;
        }
        x[i] = x0;
        y[i] = y0;
        npoints = i + 1;
        collect.push({x: x.slice(0,npoints), y: y.slice(0,npoints) })
        color.push(rgb)
        //ps_polygon(x, y, npoints, rgb, outline);
        /* second compressing part */
        i = 0;
        str = meca.NP1.strike + 180.;
        while (increment > 0. ?  str <= meca.NP1.strike + 270. + EPSIL : str >= meca.NP1.strike + 90. - EPSIL) {
            si = Math.sin(str * D2R)
            co = Math.cos(str * D2R)
            //sincosd (str, &si, &co);
            x[i] = x0 + si * radius_size;
            y[i] = y0 + co * radius_size;
            str += increment;
            i++;
        }
        x[i] = x0;
        y[i] = y0;
        npoints = i + 1;
        collect.push({x: x.slice(0,npoints), y: y.slice(0,npoints) })
        color.push(rgb)
        //ps_polygon(x, y, npoints, rgb, outline);
    }
    else {
        //console.log('other cases')
        /* other cases */
        /* first nodal plane till null axis */
        i = -1;
        increment = 1.;
        if (meca.NP1.strike > N_axis.strike)
            meca.NP1.strike -= 360.;
        str = meca.NP1.strike;
        while (Math.abs(90. - meca.NP1.dip) < EPSIL ? str <= meca.NP1.strike + EPSIL : str <= N_axis.strike + EPSIL) {
            i++;
            radius = proj_radius(meca.NP1.strike, meca.NP1.dip, str) * radius_size;
            si = Math.sin(str * D2R)
            co = Math.cos(str * D2R)
            //sincosd (str, &si, &co);
            x[i] = x0 + radius * si;
            y[i] = y0 + radius * co;
            str += increment;
        }

        /* second nodal plane from null axis */
        meca.NP2.strike += (1. + fault) * 90.;
        if (meca.NP2.strike >= 360.) meca.NP2.strike -= 360.;
        increment = fault;
        if (fault * (meca.NP2.strike - N_axis.strike) < -EPSIL) meca.NP2.strike += fault * 360.;
        str = Math.abs(90. - meca.NP2.dip) < EPSIL ? meca.NP2.strike : N_axis.strike;
        while (increment > 0. ? str <= meca.NP2.strike + EPSIL : str >= meca.NP2.strike - EPSIL) {
            i++;
            radius = proj_radius(meca.NP2.strike - (1. + fault) * 90., meca.NP2.dip, str) * radius_size;
            si = Math.sin(str * D2R)
            co = Math.cos(str * D2R)
            //sincosd (str, &si, &co);
            x[i] = x0 + radius * si;
            y[i] = y0 + radius * co;
            str += increment;
        }

        /* close the first compressing part */
        meca.NP1.strike = zero_360(meca.NP1.strike);
        meca.NP2.strike = zero_360(meca.NP2.strike);
        increment = pos_NP1_NP2 >= 0. ? -fault : fault;
        if (increment * (meca.NP1.strike - meca.NP2.strike) < - EPSIL) meca.NP1.strike += increment * 360.;
        str = meca.NP2.strike;
        while (increment > 0. ? str <= meca.NP1.strike + EPSIL : str >= meca.NP1.strike - EPSIL) {
            i++;
            si = Math.sin(str * D2R)
            co = Math.cos(str * D2R)
            //sincosd (str, &si, &co);
            x[i] = x0 + si * radius_size;
            y[i] = y0 + co * radius_size;
            str += increment;
        }

        npoints = i + 1;
        collect.push({x: x.slice(0,npoints), y: y.slice(0,npoints) })
        color.push(rgb)
        //ps_polygon(x, y, npoints, rgb, outline);

    /* first nodal plane till null axis */
        i = -1;
        meca.NP1.strike = zero_360(meca.NP1.strike + 180.);
        if (meca.NP1.strike - N_axis.strike < - EPSIL) meca.NP1.strike += 360.;
        increment = -1.;
        str = meca.NP1.strike;
        while (Math.abs(90. - meca.NP1.dip) < EPSIL ? str >= meca.NP1.strike -EPSIL : str >= N_axis.strike - EPSIL) {
            i++;
            radius = proj_radius(meca.NP1.strike - 180., meca.NP1.dip, str) * radius_size;
            si = Math.sin(str * D2R)
            co = Math.cos(str * D2R)
            //sincosd (str, &si, &co);
            x[i] = x0 + radius * si;
            y[i] = y0 + radius * co;
            str += increment;
            //console.log('str1',str, radius, meca.NP1.dip)
        }

        /* second nodal plane from null axis */
        meca.NP2.strike = zero_360(meca.NP2.strike + 180.);
        increment = -fault;
        if (fault * (N_axis.strike - meca.NP2.strike) < - EPSIL) meca.NP2.strike -= fault * 360.;
        str = Math.abs(90. - meca.NP2.dip) < EPSIL ? meca.NP2.strike : N_axis.strike;
        while (increment > 0. ? str <= meca.NP2.strike + EPSIL : str >= meca.NP2.strike - EPSIL) {
            i++;
            radius = proj_radius(meca.NP2.strike - (1. - fault) * 90., meca.NP2.dip, str) * radius_size;
            si = Math.sin(str * D2R)
            co = Math.cos(str * D2R)
            //sincosd (str, &si, &co);
            x[i] = x0 + radius * si;
            y[i] = y0 + radius * co;
            str += increment;
            //console.log('str2',str)
        }

        /* close the second compressing part */
        meca.NP1.strike = zero_360(meca.NP1.strike);
        meca.NP2.strike = zero_360(meca.NP2.strike);
        increment = pos_NP1_NP2 >= 0. ? -fault : fault;
        if (increment * (meca.NP1.strike - meca.NP2.strike) < - EPSIL) meca.NP1.strike += increment * 360.;
        str = meca.NP2.strike;
        while (increment > 0. ? str <= meca.NP1.strike + EPSIL : str >= meca.NP1.strike - EPSIL) {
            i++;
            si = Math.sin(str * D2R)
            co = Math.cos(str * D2R)
            //sincosd (str, &si, &co);
            x[i] = x0 + si * radius_size;
            y[i] = y0 + co * radius_size;
            str += increment;
            //console.log('str3',str, radius_size)
        }

        npoints = i + 1;
        collect.push({x: x.slice(0,npoints), y: y.slice(0,npoints) })
        color.push(rgb)
        //ps_polygon(x, y, npoints, rgb, outline);
    }
    return { collect, colors: color }
}

function ps_tensor(x0, y0, size,  T,  N, P, c_rgb, e_rgb, outline, plot_zerotrace) {

    const D2R = Math.PI / 180.0
    const EPSIL = 1e-5
    let b = 1;
    //GMT_LONG d, b = 1, m;
    //GMT_LONG djp, mjp;
    let j = 1
    let j2 = 0
    let j3 = 0
    lineout = 1;
    //GMT_LONG i, ii, iii, n = 0, j = 1, j2 = 0, j3 = 0;
    //GMT_LONG npoints;
    //GMT_LONG lineout = 1;
    //int rgb1[3], rgb2[3];
    //double a[3], p[3], v[3];

    let a = zeros(3)
    let p = zeros(3)
    let v = zeros(3)
    //double vi, iso, f;
    //double fir, s2alphan, alphan;
    //double cfi, sfi, can, san;
    //double cpd, spd, cpb, spb, cpm, spm;
    //double cad, sad, cab, sab, cam, sam;
    //double xz, xn, xe;
    let az = 0
    let azp = 0
    //double az = 0., azp = 0., takeoff, r;
    // double azi[3][2];
    // double x[400], y[400], x2[400], y2[400], x3[400], y3[400];
    // double xp1[800], yp1[800], xp2[400], yp2[400];
    // double radius_size;
    // double si, co;
    // int jp_flag;
    // int bigisotestv0, bigisotestv2;

    let collect = []
    let colors = []
    let azi = zeros([3,2])

    let n = 0

    x = zeros(400)
    y = zeros(400)
    x2 = zeros(400)
    y2 = zeros(400)
    x3 = zeros(400)
    y3 = zeros(400)

    xp1 = zeros(800)
    yp1 = zeros(800)
    xp2 = zeros(400)
    yp2 = zeros(400)

    //console.log(T,N,P);
    a[0] = T.strike; a[1] = N.strike; a[2] = P.strike;
    p[0] = T.dip; p[1] = N.dip; p[2] = P.dip;
    v[0] = T.val; v[1] = N.val; v[2] = P.val;

    vi = (v[0] + v[1] + v[2]) / 3.;
    for (i=0; i<=2; i++)
        v[i] = v[i] - vi;

    /*Redundant cases for pure explosion and large isotropic with no nodes removed*/
    /*DS Dreger 2011*/

    /*Test to determine if angle exceeds plotting capability*/
    /*Used to choose the dominant eigenvalue after Frohlich for plotting purposes*/
    /*DS Dreger 2011*/
    bigisotestv0=0;
    bigisotestv2=0;
    for (i=0; i<360; i++) {
        fir =  i * D2R;
        f = -v[1]/v[0];
        iso = vi/v[0];
        s2alphan = (2. + 2. * iso) / (3. + (1. - 2. * f) * Math.cos(2. * fir));
        if(s2alphan > 1.0)
            bigisotestv0 += 1;
        f = -v[1]/v[2];
        iso = vi/v[2];
        s2alphan = (2. + 2. * iso) / (3. + (1. - 2. * f) * Math.cos(2. * fir));
        if(s2alphan > 1.0)
            bigisotestv2 += 1;
    }
    /*End Test*/

    let radius_size = size * 0.5;

    /*Determine the eigenvalue to use based on s2alphan test above*/
    /*Define the shading convention*/
    /*DS Dreger 2011*/
    if(bigisotestv0 == 0) {
        d=0;
        m=2;
        rgb1 = c_rgb
        rgb2 = e_rgb
    } else if(bigisotestv2 == 0) {
        d=2;
        m=0;
        rgb1 = e_rgb
        rgb2 = c_rgb
    } else {
        //fprintf(stderr,"Error bigisotest failed, exiting\n");
        console.log("Error bigisotest failed")
        return undefined
    }

    if (plot_zerotrace)
        vi = 0.;

    f = - v[1] / v[d];
    iso = vi / v[d];
    jp_flag = 0;
    djp = -1;
    mjp = -1;


    //sincosd (p[d], &spd, &cpd);
    //sincosd (p[b], &spb, &cpb);
    //sincosd (p[m], &spm, &cpm);
    let spd = Math.sin(p[d] * D2R)
    let cpd = Math.cos(p[d] * D2R)
    let spb = Math.sin(p[b] * D2R)
    let cpb = Math.cos(p[b] * D2R)
    let spm = Math.sin(p[m] * D2R)
    let cpm = Math.cos(p[m] * D2R)

    //sincosd (a[d], &sad, &cad);
    //sincosd (a[b], &sab, &cab);
    //sincosd (a[m], &sam, &cam);
    let sad = Math.sin(a[d] * D2R)
    let cad = Math.cos(a[d] * D2R)
    let sab = Math.sin(a[b] * D2R)
    let cab = Math.cos(a[b] * D2R)
    let sam = Math.sin(a[m] * D2R)
    let cam = Math.cos(a[m] * D2R)

    //console.log("f",f,iso)
    //console.log("p", spd, spb, spm, cpd, cpb, cpm);
    //console.log("a", sad, sab, sam, cad, cab, cam);

    for (i=0; i<360; i++) {
        fir = i * D2R;
        s2alphan = (2. + 2. * iso) / (3. + (1. - 2. * f) * Math.cos(2. * fir));

        alphan = Math.asin(Math.sqrt(s2alphan));
        let sfi = Math.sin(fir)
        let cfi = Math.cos(fir)
        //sincos (fir, &sfi, &cfi);
        let san = Math.sin(alphan)
        let can = Math.cos(alphan)
        //sincos (alphan, &san, &can);

        xz = can * spd + san * sfi * spb + san * cfi * spm;
        xn = can * cpd * cad + san * sfi * cpb * cab + san * cfi * cpm * cam;
        xe = can * cpd * sad + san * sfi * cpb * sab + san * cfi * cpm * sam;

        if (Math.abs(xn) < EPSIL && Math.abs(xe) < EPSIL) {
            takeoff = 0.;
            az = 0.;
        } else {
            az = Math.atan2(xe, xn);
            if (az < 0.)
                az += Math.PI * 2.;
            takeoff = Math.acos(xz / Math.sqrt(xz * xz + xn * xn + xe * xe));
        }
        if (takeoff > (Math.PI / 2.0)) {
            takeoff = Math.PI - takeoff;
            az += Math.PI;
            if (az > Math.PI * 2.)
                az -= Math.PI * 2.;
        }
        r = Math.sqrt(2.0) * Math.sin(takeoff / 2.);
        let si = Math.sin(az)
        let co = Math.cos(az)
        //sincos (az, &si, &co);
        if (i == 0) {
            azi[i][0] = az;
            x[i] = x0 + radius_size * r * si;
            y[i] = y0 + radius_size * r * co;
            azp = az;
        } else {
            //console.log('az',i, azp, az, az-azp, D2R*10)
            if (Math.abs(Math.abs(az - azp) - Math.PI) < D2R * 10. && takeoff > 80 * D2R) {
                //console.log('*******n',i, n,n+1, az, azp, Math.abs(Math.abs(az - azp) - Math.PI), 10.0 * D2R)
                azi[n][1] = azp;
                n += 1
                azi[n][0] = az;
            }
            if (Math.abs(Math.abs(az -azp) - Math.PI * 2.) < D2R * 2.) {
                if (azp < az) {
                    azi[n][0] += Math.PI * 2.;
                } else {
                    azi[n][0] -= Math.PI * 2.;
                }
            }
            if(n == 0) {
                x[j] = x0 + radius_size * r * si;
                y[j] = y0 + radius_size * r * co;
                j++;
            } else if (n == 1) {
                x2[j2] = x0 + radius_size * r * si;
                y2[j2] = y0 + radius_size * r * co;
                j2++;
            } else if (n == 2) {
                x3[j3] = x0 + radius_size * r * si;
                y3[j3] = y0 + radius_size * r * co;
                j3++;
            }
        }
        azp = az;
    }
    azi[n][1] = az;

    //console.log('n',n)
    //console.log('azi',azi)
    //console.log('color', rgb2)
    collect.push( xcircle([x0, y0], radius_size ))
    colors.push(rgb2)
    //ps_circle(x0, y0, radius_size*2., rgb2, lineout);

    if(n == 0) {
        let i
        //for (i=0; i<360; i++) {
        //    xp1[i] = x[i];
        //    yp1[i] = y[i];
        //}
        npoints = 360;
        //ps_polygon(xp1, yp1, npoints, rgb1, outline);
        collect.push({ x: x.slice(0, npoints), y: y.slice(0, npoints) })
        colors.push(rgb1)
    }
    if(n == 1) {
        //case 1 :
        for (i=0; i<j; i++) {
            xp1[i] = x[i]; yp1[i] = y[i];
        }
        if (azi[0][0] - azi[0][1] > Math.PI)
            azi[0][0] -= Math.PI * 2.;
        else if (azi[0][1] - azi[0][0] > Math.PI)
            azi[0][0] += Math.PI * 2.;

        if (azi[0][0] < azi[0][1]) {
            for (az = azi[0][1] - D2R; az > azi[0][0]; az -= D2R) {
                //sincos (az, &si, &co);
                let si = Math.sin(az);
                let co = Math.cos(az);
                xp1[i] = x0 + radius_size * si;
                yp1[i] = y0 + radius_size * co;
                i += 1
            }
        } else {
            for (az = azi[0][1] + D2R; az < azi[0][0]; az += D2R) {
                //sincos (az, &si, &co);
                let si = Math.sin(az);
                let co = Math.cos(az);
                xp1[i] = x0 + radius_size * si;
                yp1[i] = y0 + radius_size * co;
                i += 1
            }
        }
        npoints = i;
        //ps_polygon(xp1, yp1, npoints, rgb1, outline);
        collect.push({ x: xp1.slice(0, npoints), y: yp1.slice(0, npoints) })
        colors.push(rgb1)

        for (i=0; i<j2; i++) {
            xp2[i] = x2[i]; yp2[i] = y2[i];
        }
        if (azi[1][0] - azi[1][1] > Math.PI)
            azi[1][0] -= Math.PI * 2.;
        else if (azi[1][1] - azi[1][0] > Math.PI)
            azi[1][0] += Math.PI * 2.;

        if (azi[1][0] < azi[1][1]) {
            for (az = azi[1][1] - D2R; az > azi[1][0]; az -= D2R) {
                //sincos (az, &si, &co);
                let si = Math.sin(az);
                let co = Math.cos(az);
                xp2[i] = x0 + radius_size * si;
                yp2[i] = y0 + radius_size * co;
                i += 1
            }
        } else {
            for (az = azi[1][1] + D2R; az < azi[1][0]; az += D2R) {
                //sincos (az, &si, &co);
                let si = Math.sin(az);
                let co = Math.cos(az);
                xp2[i] = x0 + radius_size * si;
                yp2[i] = y0 + radius_size * co;
                i += 1
            }
        }
        npoints = i;
        //ps_polygon(xp2, yp2, npoints, rgb1, outline);
        collect.push({ x: xp2.slice(0, npoints), y: yp2.slice(0, npoints) })
        colors.push(rgb1)
    }
    if(n == 2) {
        //case 2 :
        let i
        for (i=0; i<j3; i++) {
            xp1[i] = x3[i];
            yp1[i] = y3[i];
        }
        //console.log("0 ", xp1[0],yp1[0])
        //console.log("1 ", xp1[1],yp1[1])
        //console.log(i-3, xp1[i-3],yp1[i-3])
        //console.log(i-2, xp1[i-2],yp1[i-2])
        //console.log(i-1, xp1[i-1],yp1[i-1])
        for (ii=0; ii<j; ii++) {
            xp1[i] = x[ii];
            yp1[i] = y[ii];
            i += 1
        }
        //console.log(i-2, xp1[i-2],yp1[i-2])
        //console.log(i-1, xp1[i-1],yp1[i-1])

        if (azi[2][0] - azi[0][1] > Math.PI)
            azi[2][0] -= Math.PI * 2.;
        else if (azi[0][1] - azi[2][0] > Math.PI)
            azi[2][0] += Math.PI * 2.;
        if (azi[2][0] < azi[0][1]) {
            //console.log(`AZ ${azi[2][0].toFixed(4)} < ${azi[0][1].toFixed(4)}`)
            for (az = azi[0][1] - D2R; az > azi[2][0]; az -= D2R) {
                //sincos (az, &si, &co);
                let si = Math.sin(az);
                let co = Math.cos(az);
                xp1[i] = x0+ radius_size * si;
                yp1[i] = y0+ radius_size * co;
                i += 1
            }
        } else {
            for (az = azi[0][1] + D2R; az < azi[2][0]; az += D2R) {
                //sincos (az, &si, &co);
                let si = Math.sin(az);
                let co = Math.cos(az);
                xp1[i] = x0+ radius_size * si;
                yp1[i] = y0+ radius_size * co;
                i += 1
            }
        }
        //console.log(i-2, xp1[i-2],yp1[i-2])
        //console.log(i-1, xp1[i-1],yp1[i-1])

        npoints = i;
        //for(let i = 0; i < npoints; i++) {
        //    console.log(i, xp1[i].toFixed(4), yp1[i].toFixed(4))
        //}
        collect.push({ x: xp1.slice(0, npoints), y: yp1.slice(0, npoints) })
        colors.push(rgb1)
        //ps_polygon(xp1, yp1, npoints, rgb1, outline);
        for (i=0; i<j2; i++) {
            xp2[i] = x2[i]; yp2[i] = y2[i];
        }
        if (azi[1][0] - azi[1][1] > Math.PI)
            azi[1][0] -= Math.PI * 2.;
        else if (azi[1][1] - azi[1][0] > Math.PI)
            azi[1][0] += Math.PI * 2.;

        if (azi[1][0] < azi[1][1]) {
            for (az = azi[1][1] - D2R; az > azi[1][0]; az -= D2R) {
                //sincos (az, &si, &co);
                let si = Math.sin(az);
                let co = Math.cos(az);
                xp2[i] = x0+ radius_size * si;
                yp2[i] = y0+ radius_size * co;
                i += 1
            }
        } else {
            for (az = azi[1][1] + D2R; az < azi[1][0]; az += D2R) {
                //sincos (az, &si, &co);
                let si = Math.sin(az);
                let co = Math.cos(az);
                xp2[i] = x0+ radius_size * si;
                yp2[i] = y0+ radius_size * co;
                i += 1
            }
        }
        npoints = i;
        //ps_polygon(xp2, yp2, npoints, rgb1, outline);
        collect.push({ x: xp2.slice(0, npoints), y: yp2.slice(0, npoints)})
        colors.push(rgb1)
        //break;
    }
    //return(radius_size*2.);
    return { collect, colors }
}


L.BeachballIcon = L.Icon.extend({
    options: {

    },

    collect: [],
    colors: [],
    size: 200,

    createIcon: function() {
        var e = document.createElement("canvas");
        this._setIconStyles(e, "icon");
        var size = this.options.iconSize;
        const out = plotmt(this.options.mech, {
            size: size[0],
            tcolor: this.options.tcolor,
            pcolor: this.options.pcolor
        })
        this.collect = out.collect;
        this.colors = out.colors;
        e.width = size[0];
        e.height = size[1];
        this.ctx = e.getContext("2d");
        this.draw(this.ctx, size[0], size[1]);

        return e;
    },

    draw: function(ctx, w, h) {
        if(!ctx) return;
        ctx.clearRect(0, 0, w, h);
        ctx.save()
        ctx.translate(w/2, h/2);
        ctx.beginPath();
        for(let j = 0; j < this.collect.length; j++) {
            x = this.collect[j].x
            y = this.collect[j].y
            ctx.beginPath()
            ctx.fillStyle = this.colors[j]
            ctx.strokeStyle = this.colors[j]//'rgba(0,0,0,1.0)'
            ctx.lineWidth = 0
            ctx.moveTo(x[0], -y[0] )
            for(let i = 1; i < x.length; i++) {
                ctx.lineTo(x[i], -y[i])
            }
            ctx.closePath()
            ctx.fill()
            ctx.stroke()
        }
    },
});

L.BeachballMarker = L.Marker.extend({
    options: {},
    initialize: function(pos, options) {
        this.setLatLng(pos)
        L.setOptions(this, options)
    },
    pos: function () {
        return this.getLatLng();
    },
});

L.beachballMarker = function(pos, options) {
    const iconSizeDefault = 20;
    options.tcolor = ("tcolor" in options) ? options.tcolor : 'black'
    options.pcolor = ("pcolor" in options) ? options.pcolor : 'white'
    options.iconSize = ("iconSize" in options) ? options.iconSize : [iconSizeDefault,iconSizeDefault]
    options.mech = ("mech" in options) ? options.mech : {xx:0,yy:0,zz:0,xy:0,xz:0,yz:1}

    options.iconSize = ("size" in options) ? options.size : [iconSizeDefault,iconSizeDefault]
    options.iconAnchor = ("iconAnchor" in options) ? options.iconAnchor : [options.iconSize[0]/2, options.iconSize[1]/2]

    options.icon = new L.BeachballIcon(options);

    return new L.BeachballMarker(pos, options)
};
