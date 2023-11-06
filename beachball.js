"use strict;"

let DEBUG = false;
const log = function() {
    if(DEBUG)
        console.log.apply(console, arguments)
    return arguments[0]
}
export function set_debug(val) { DEBUG = val; }

function zeros(n) {
    if(Array.isArray(n) && n.length == 2) {
        let out = []
        for(let i = 0; i < n[0]; i++) {
            let row = []
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

function xcircle(center, radius, inc = 10) {
    let x = []
    let y = []
    for(let i = 0; i <= 360; i += inc) {
        x.push(center[0] + radius * Math.cos(i * Math.PI/180.0))
        y.push(center[1] + radius * Math.sin(i * Math.PI/180.0))
    }
    return {x, y}
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
    console.log("Too many iterations in routine jacobi");
    return undefined
}

function ROTATE(a,i,j,k,l, s, tau) {
    let g = a[i-1][j-1]
    let h = a[k-1][l-1]
    a[i-1][j-1] = g - s * (h + g * tau)
    a[k-1][l-1] = h + s * (g - h * tau)
}

class Pole {
    n = [0,0,0]
    constructor(n) {
        const R2D = 180.0 / Math.PI
        this.n = norm(n)

        // Convert to Angles
        const xy = Math.sqrt(n[0]*n[0] + n[1]*n[1])
        // Really trend and plunge
        const flip = (this.n[2] > 0) ? -1 : 1
        const m = this.n.map(v => v * flip)

        this.dip = Math.asin(-m[2]) * R2D
        this.strike = 90 - Math.atan2(m[1], m[0]) * R2D
    }
    plane() {
        const p = normal_plane(this.n)
        return new Plane(p.str, p.dip, this.n)
    }
}

function sign(val) {
    return (val < 0) ? -1 : 1
}

export class Plane {
    strike = 0
    dip = 90
    constructor(strike, dip, normal = undefined) {
        this.strike = zero_360(strike)
        this.dip = dip
        this.n = (normal) || this.normal().n
    }
    normal() {
        return new Pole(plane_normal(this.strike, this.dip))
    }
    strike_vec() {
        return normal_to_strike_vec(this.normal)
    }
    is_horizontal() {
        return close(this.n[2], 1) || close(this.n[2], -1)
    }
    aux_plane(rake) {
        if(this.is_horizontal()) {
            // Horizontal Planes have vertical auxillary planes
            return new Plane(-sign(rake) * 90 + this.strike - rake, 90)
        }
        const n2 = new Pole(normal_second_plane(this.n, rake, this.strike))
        return n2.plane()
    }
    to_nodal_plane(rake) {
        return new NodalPlane(this.strike, this.dip, rake, this.n)
    }
}
export class NodalPlane extends Plane {
    rake = 0
    constructor(strike, dip, rake, normal = undefined) {
        super(strike, dip, normal)
        this.rake = rake
        this.rake_vector = normal_second_plane(this.n, rake, this.strike)
    }
    aux_plane() {
        let p = super.aux_plane(this.rake)
        if(p.is_horizontal()) {
            p = new Plane(this.strike + 180, 0)
        }
        const ax = this.stress_axes(p)
        const r = rake3(this.n, p.n, ax.t.n, ax.p.n, this.strike, p.strike)
        return new NodalPlane(p.strike, p.dip, r.rake2, p.n)
    }
    stress_axes(np) {
        const b = norm([ ... np.n ])
        let t = [this.n[0] + b[0], this.n[1] + b[1], this.n[2] + b[2]]
        let p = [this.n[0] - b[0], this.n[1] - b[1], this.n[2] - b[2]]

        let sa = stress_axis_hanging_wall(this.n, t, p)
        let poss = norm(vsub(sa.t, sa.p))
        const d = dot(poss, this.rake_vector);
        if(d < 0) { // Rake is opposite from prescribed rake, flip P and T Axis
            [t,p] = [p,t]
        }
        t = new Pole(t)
        p = new Pole(p)
        return {t, p}
    }
}

function vadd(a,b) {
    return [0,1,2].map(i => a[i] + b[i])
}
function vsub(a,b) {
    return [0,1,2].map(i => a[i] - b[i])
}

function computeM0(M) {
    let sum = 0.0
    for(let i = 0 ; i < 3; i++) {
        for(let j = 0; j < 3; j++) {
            sum += M[i][j]
        }
    }
    return 1.0 / Math.sqrt(2.0) * Math.sqrt(sum)
}

export function mt2axes(mt) {
    const EPSIL = 1e-5
    /*
    Calculates the principal axes of a given moment tensor.

    :param mt: :class:`~MomentTensor`
    :return: tuple of :class:`~PrincipalAxis` T, N and P

    Adapted from ps_tensor / utilmeca.c /
    `Generic Mapping Tools (GMT) <https://gmt.soest.hawaii.edu>`_.
    */
    const np = 3
    const R2D = 180.0 / Math.PI;
    let M = []
    if('xx' in mt) {
        M = [ [mt.xx, mt.xy, mt.xz],
              [mt.xy, mt.yy, mt.yz],
              [mt.xz, mt.yz, mt.zz] ]
    }
    if('rr' in mt) {
        M = [ [mt.rr, mt.rt, mt.rp],
              [mt.rt, mt.tt, mt.tp],
              [mt.rp, mt.tp, mt.pp] ]
    }
    const M0 = computeM0(M);
        //(d, v) = np.linalg.eigh(M)
    //console.log('M',JSON.stringify(M))
    //console.log('M',M)
    let ret = jacobi23(M)
    if(ret == undefined) {
        console.error('Error in mt2axes', M)
        return undefined;
    }
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
    let t = { val: d[0], strike: az[0], dip: pl[0], v: [-v[2][0], v[1][0], -v[0][0]] }
    let n = { val: d[1], strike: az[1], dip: pl[1], v: [-v[2][1], v[1][1], -v[0][1]] }
    let p = { val: d[2], strike: az[2], dip: pl[2], v: [-v[2][2], v[1][2], -v[0][2]] }
    //console.log('t',t)
    //console.log('n',n)
    //console.log('p',p)
    return {t, n, p}
}

function plane_from_pole(n) {
    const R2D = 180.0 / Math.PI
    if(n[2] < 0)
        n = [-n[0],-n[1],-n[2]]
    // Angle of the Normal vector from the vertical = Angle of the Plane from horizontal
    let dip = R2D * Math.atan2(Math.sqrt(n[0]*n[0] + n[1]*n[1]), n[2])
    let strike = zero_360(R2D * Math.atan2(n[1], -n[0]))
    return { strike, dip }
}

function len(n) {
    return Math.sqrt(dot(n,n))
}

export function norm(v) {
    const vlen = len(v)
    return v.map(x => x / vlen)
}

function dot(a,b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
}

function between(start0, mid0, end0) {
    //console.log('between', start0, mid0, end0)
    let end = zero_360(end0) - zero_360(start0)
    let mid = zero_360(mid0) - zero_360(start0)
    //console.log('       ', end, mid)
    end = (end < 0) ? end + 360 : end
    mid = (mid < 0) ? mid + 360 : mid
    //console.log('       ', end, mid)
    return mid <= end
}

function close(a, b, thresh = 1e-5) {
    return Math.abs(a-b) < thresh
}

function stress_axis_hanging_wall(n, t, p) {
    // Get complementary P and T axes
    const t1 = norm([...t])
    const t2 = norm(t.map(v => -v))
    const p1 = norm([...p])
    const p2 = norm(p.map(v => -v))

    let n1_t1 = dot(n,t1)
    let n1_p1 = dot(n,p1)

    let t0 = (n1_t1 < 0) ? t1 : t2;
    let p0 = (n1_p1 < 0) ? p1 : p2;

    return {t: t0, p: p0}

}

function rake3(n1, n2, t, p, str1, str2) {
    const EPSIL = 1e-5
    let R2D = 180.0 / Math.PI
    let D2R = Math.PI / 180.0
    let s1 = normal_to_strike_vec(n1)
    let s2 = normal_to_strike_vec(n2)
    //console.log('RAKE3',n1,n2)
    // Strike should be from downward facing normal

    if(close(n1[2],1) || close(n1[2],-1))
        s1 = norm([Math.cos((90-str1) * D2R), Math.sin((90-str1) * D2R), 0])
    if(close(n2[2],1) || close(n2[2],-1))
        s2 = norm([Math.cos((90-str2) * D2R), Math.sin((90-str2) * D2R), 0])
    if(n1[2] > EPSIL) {
        s1 = norm([-s1[0], -s1[1], 0])
    }
    if(n2[2] > EPSIL) {
        s2 = norm([-s2[0], -s2[1], 0])
    }

    // Get complementary P and T axes
    const t1 = norm([...t])
    const t2 = norm(t.map(v => -v))
    const p1 = norm([...p])
    const p2 = norm(p.map(v => -v))

    log("T", t1)
    log("P", p1)

    // Compute angle between plane normals and P and T axes
    //   We want to chose the P and T axes on the side opposite from
    //   the plane normal (this is the hanging wall side)
    let n1_t1 = dot(n1,t1)
    let n1_p1 = dot(n1,p1)
    let n2_t1 = dot(n2,t1)
    let n2_p1 = dot(n2,p1)

    let T1 = (n1_t1 < 0) ? t1 : t2
    let P1 = (n1_p1 < 0) ? p1 : p2
    let T2 = (n2_t1 < 0) ? t1 : t2
    let P2 = (n2_p1 < 0) ? p1 : p2

    // Compute the slip direction on the plane from the P and T axes
    //   Assume the P axis is pointing in and T axis is pointing out
    //   Direction should be either along the other plane's normal or in
    //   the opposite direction
    //   Either way we can compute the slip from this direction.
    let R1 = norm([0,1,2].map(i => -P1[i] + T1[i]))
    let R2 = norm([0,1,2].map(i => -P2[i] + T2[i]))

    log('R1',_fa(R1), _fa(n2), dot(R1,n2))
    log('R2',_fa(R2), _fa(n1), dot(R2,n1))

    let arg1 = dot(R1,s1) / (len(R1)*len(s1))
    let arg2 = dot(R2,s2) / (len(R2)*len(s2))
    arg1 = Math.max(Math.min(arg1, 1.0), -1.0)
    arg2 = Math.max(Math.min(arg2, 1.0), -1.0)

    log('N,S1', _fa(n1), _fa(s1), zero_360(90- R2D * Math.atan2(s1[1], s1[0])))
    log('N,S2', _fa(n2), _fa(s2), zero_360(90- R2D * Math.atan2(s2[1], s2[0])))
    log()
    log('D N2 S1', _fa(R1), _fa(s1))
    log('D N2,S1',  dot(R1,s1), (len(n2) * len(s1)))
    log('D N2,S1',  arg1)
    log('D N2,S1', R2D * Math.acos( arg1) )
    log()
    log('D N1 S2', _fa(R2), _fa(s2))
    log('D N1,S2',  dot(R2,s2))
    log('D N1,S2',  arg2 )
    log('D N1,S2', R2D * Math.acos( arg2 ))

    let rake1 = R2D * Math.acos(arg1)
    let rake2 = R2D * Math.acos(arg2)

    // -180 -   0 - 180  (Normal Strike)
    //    0 - 180 - 360  (Flipped Strike)
    // Get Correct Orientation
    if(R1[2] < -EPSIL)
        rake1 *= -1
    if(R2[2] < -EPSIL)  {// Downward pointing normal should be a negative angle
        rake2 *= -1   // Clockwise is negative
        log("Downward pointing normal n1, flip rake sign on n1", rake2)
    }
    
    if(n1[2] > EPSIL) // Upward pointing normal will have a flipped strike
        rake1 += 180
    if(n2[2] > EPSIL) {
        rake2 += 180
        log("Upward   pointing normal n2, flip rake by 180 on n2", rake2)
    }

    // Horizontal Planes
    if(close(n2[2],1) || close(n2[2],-1) ) {
        //const rp = rake2 * D2R
        //const r2p = [s2[0] * Math.cos(rp) - s2[1] * Math.sin(rp),
        //             s2[0] * Math.sin(rp) + s2[1] * Math.cos(rp),
        //             0]
        // Horizontal rotation from Strike to Rake
        const rm = -rake2 * D2R
        const r2m = [s2[0] * Math.cos(rm) - s2[1] * Math.sin(rm),
                     s2[0] * Math.sin(rm) + s2[1] * Math.cos(rm),
                     0]
        if(close(dot(R2,r2m),1))
            rake2 *= -1
    }
    if(close(n1[2],1) || close(n1[2],-1)) {
        const rm = -rake1 * D2R
        const r1m = [s1[0] * Math.cos(rm) - s1[1] * Math.sin(rm),
                     s1[0] * Math.sin(rm) + s1[1] * Math.cos(rm),
                     0]
        if(close(dot(R1,r1m),1))
            rake1 *= -1
    }

    rake1 = z180(rake1)
    rake2 = z180(rake2)
    //console.log(` RAKE1 ${f(rake1)} ${R1[2] > EPSIL}`)
    //console.log(` RAKE2 ${f(rake2)} ${R2[2] > EPSIL} `)

    return {rake1, rake2}
}

export function plotmt(M, options = {}) {
    const size = options.size || 100
    const tcolor = options.tcolor || 'black'
    const pcolor = options.pcolor || 'white'
    const EPSILON = 1e-5

    const ax = mt2axes(M)
    if(close(ax.t.val,0) && close(ax.p.val, 0) && close(ax.n.val, 0))
        return undefined;
    if(ax == undefined) {
        console.error('Error in plotmt',M)
        return undefined;
    }
    //console.log('AX',ax)
    if(Math.abs(ax.n.val) < EPSILON && Math.abs(ax.t.val + ax.p.val) < EPSILON ) {
        let mp = axes_to_planes_dc(ax.t, ax.p)
        mp.ax = new NodalPlane(mp.NP1.strike, mp.NP1.dip, mp.NP1.rake).stress_axes(mp.NP2)
        let out =  ps_mechanism(0,0,structuredClone(mp),size, tcolor,pcolor,1)
        out.ax = mp.ax
        out.dc = true
        out.iso = false
        out.clvd = false
        out.size = size
        out.np = mp
        return out
    } else {
        log(M.rr,M.tt,M.pp,M.rt,M.rp, M.tp)
        //console.log("tensor", ax.n.val, ax.t.val + ax.p.val)
        let out = ps_tensor(0,0, size, ax.t, ax.n, ax.p, tcolor, pcolor, true, false );
        if(!out) {
            console.log("out undefined", M, ax)
            return undefined
        }
        out.ax = ax
        out.dc = false
        out.iso = Math.abs(ax.t.val - ax.p.val) < EPSILON &&
            Math.abs(ax.t.val - ax.n.val) < EPSILON
        out.clvd = is_clvd(ax)
        out.n_zero = close(ax.n.val,0) && !close(ax.t.val,0) && !close(ax.p.val,0)
        out.np_zero = close(ax.n.val,0) && !close(ax.t.val,0) && close(ax.p.val,0)
        out.nt_zero = close(ax.n.val,0) && close(ax.t.val,0) && !close(ax.p.val,0)
        out.size = size
        return out
    }
}

function is_clvd(ax) {
    const EPSIL = 1e-5
    let t = ax.t.val
    let p = ax.p.val
    let n = ax.n.val
    if(close(t + 2*p,0) && close(t + 2*n, 0))
        return true;
    if(close(p + 2*t, 0) && close(p + 2*n, 0))
        return true;
    const avg = (t + p + n) / 3
    t -= avg
    n -= avg
    p -= avg
    if(close(t + 2*p, 0) && close(t + 2*n, 0))
        return true;
    if(close(p + 2*t, 0) && close(p + 2*n, 0))
        return true;
    return false
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
    //console.log('null_axis_dip',den, dip1, dip2)
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

    let sd1 = Math.sin(dip1 * D2R)
    let cd1 = Math.cos(dip1 * D2R)
    let sd2 = Math.sin(dip2 * D2R)
    let cd2 = Math.cos(dip2 * D2R)

    let ss1 = Math.sin(str1 * D2R)
    let cs1 = Math.cos(str1 * D2R)
    let ss2 = Math.sin(str2 * D2R)
    let cs2 = Math.cos(str2 * D2R)

    let cosphn = sd1 * cs1 * cd2 - sd2 * cs2 * cd1;
    let sinphn = sd1 * ss1 * cd2 - sd2 * ss2 * cd1;
    if (Math.sin((str1 - str2) * D2R) < 0.) {
        cosphn = -cosphn;
        sinphn = -sinphn;
    }
    let phn = Math.atan2(sinphn, cosphn) * R2D
    if (phn < 0.0)
        phn += 360.;
    return phn
}

function proj_radius(str1, dip1, str) {
    /*
      Compute the vector radius for a given strike,
      equal area projection, inferior sphere.
      Strike and dip of the plane are given.

      Genevieve Patau
    */

    const D2R = Math.PI / 180.0
    const EPSIL = 1e-5
    let dip = 0
    let r = 1;
    if (Math.abs(dip1 - 90.) < EPSIL) {
      r = (Math.abs(str - str1) < EPSIL || Math.abs(str - str1 - 180) < EPSIL) ? 1. : 0.;
    } else {
        dip = Math.atan(Math.tan(dip1 * D2R) * Math.sin((str - str1) * D2R));
        r = Math.sqrt(2.) * Math.sin(Math.PI/4.0 - dip/2.);
    }
    return r;
}

export function z180(rake) {
    if(rake > 180)
        rake -= 360
    if(rake <= -180.0)
        rake += 360.0
    return rake
}
export function zero_360(str) {
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
    let sp2 = -am * temp
    temp = ss * cr
    temp -= sr * cs * cd1
    let cp2 = am * temp
    let str2 = Math.atan2(sp2, cp2) * R2D;
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
    return Math.atan2(sinrake2, -am * sd * ss) * R2D
}

export function define_second_plane(np1) {
    let strike = zero_360(computed_strike1(np1))
    let dip = computed_dip1(np1)
    let rake = computed_rake1(np1, strike, dip)
    return { strike, dip, rake }
}


function strike_slip_taxis(meca) {
    const s1 = meca.NP1.strike
    const t = meca.ax.t.strike
    if(between(s1, t, s1 + 90)) {
        return 1
    }
    if(between(s1 + 180, t, s1 + 270)) {
        return 1
    }
    if(between(s1 - 90, t, s1)) {
        return -1
    }
    if(between(s1 + 90, t, s1 + 180)) {
        return -1
    }
    console.error("Error finding taxis for strike slip fault", s1, t)
    return 1
}

function arc(x0, y0, radius_size, str0, str1) {
    const D2R = Math.PI / 180.0
    let x = []
    let y = []
    let str = str0
    let increment = (str0 < str1) ? 1 : -1
    while((increment > 0) ? str <= str1 : str >= str1) {
        let si = Math.sin(str * D2R)
        let co = Math.cos(str * D2R)
        x.push( x0 + si * radius_size )
        y.push( x0 + co * radius_size )
        str += increment
    }
    return {x, y}
}

export function ps_mechanism(x0, y0, meca, size, rgb, ergb, outline)  {

    const D2R = Math.PI / 180.0
    const EPSIL = 1e-5

    let collect = []
    let color = []
    console.log("PS MECHANISM")
    let pos_NP1_NP2 = Math.sin((meca.NP1.strike-meca.NP2.strike)*D2R)

    if(meca.ax === undefined) {
        meca.ax = dc2axe(meca)
    }

    let fault = 1
    //console.log('NP1', meca.NP1)
    //console.log('NP2', meca.NP2)
    // Fault +1 - Compressional, -1 Extensional
    //console.log('MECA', meca)
    //if(Math.abs(meca.NP1.rake) < EPSIL)
    if(close(meca.NP1.rake,0) || close(meca.NP1.rake,180) || close(meca.NP1.rake,-180))
        fault = meca.NP2.rake / Math.abs(meca.NP2.rake)
    else
        fault = meca.NP1.rake / Math.abs(meca.NP1.rake)
    //console.log('fault', fault, 'rake1,2', _f(meca.NP1.rake), _f(meca.NP2.rake), _f(meca.NP1.strike), _f(meca.NP2.strike))
    //struct AXIS N_axis;
    let N_axis = { strike: 0, dip: 0, rake: 0};
    /* compute null axis strike and dip */
    N_axis.dip = null_axis_dip(meca.NP1.strike, meca.NP1.dip, meca.NP2.strike, meca.NP2.dip);
    if (Math.abs(90. - N_axis.dip) < EPSIL)
        N_axis.strike = meca.NP1.strike;
    else
        N_axis.strike = null_axis_strike(meca.NP1.strike, meca.NP1.dip, meca.NP2.strike, meca.NP2.dip);

    /* compute radius size of the bubble */
    let radius_size = size * 0.5;

    log("null axis", N_axis.strike, N_axis.dip)
    /* outline the bubble */
    //ps_plot(x0 + radius_size, y0, 3);
    /*  argument is DIAMETER!!*/
    //ps_circle(x0, y0, radius_size*2., ergb, lineout);
    let x = []
    let y = []
    collect.push( xcircle([x0,y0], radius_size) )
    color.push(ergb)

    if(close(meca.NP1.dip,0) && close(meca.NP2.dip, 90) ||
       close(meca.NP2.dip,0) && close(meca.NP1.dip, 90)) {
        console.log("HORIZONTAL AND VERTICAL PLANES")
        let str = (close(meca.NP1.dip,90)) ? meca.NP1.strike : meca.NP2.strike
        let increment = between(str, meca.ax.t.strike, str + 180) ? 1 : -1
        let i = 0
        console.log('STRIKE', str, increment, str+180*increment)
        let pts = arc(x0, y0, radius_size, str, str + 180 * increment)
        pts.x.push( pts.x[0] )
        pts.y.push( pts.y[0] )
        console.log(pts)
        collect.push( { x: pts.x, y: pts.y })
        color.push(rgb)
        return {collect, colors: color }
    } else if (Math.abs(pos_NP1_NP2) < EPSIL) {
        console.log("pure compression/extensional fault")
        /* pure normal or inverse fault (null axis strike is determined
           with + or - 180 degrees. */
        /* first nodal plane part */
        let t_axis_vertical = meca.ax.t.dip > meca.ax.p.dip
        if(close(meca.ax.t.dip, meca.ax.p.dip)) { // Vertical Dip Slip
            t_axis_vertical = (close(meca.NP1.rake,0)) ? meca.NP2.rake > 0 : meca.NP1.rake > 0;
            console.log('DIPS EQUAL')
        }
        console.log('MECA',meca)
        console.log('t_axis_vertical', t_axis_vertical)
        let i = -1;
        let increment = 1.;
        let str = meca.NP1.strike;
        while (str <= meca.NP1.strike + 180. + EPSIL) {
            i++;
            let radius = proj_radius(meca.NP1.strike, meca.NP1.dip, str) * radius_size;
            //sincosd (str, &si, &co);
            let si = Math.sin(str * D2R)
            let co = Math.cos(str * D2R)
            x[i] = x0 + radius * si;
            y[i] = y0 + radius * co;
            str += increment;
        }
        if (!t_axis_vertical) {
            /* normal fault, close first compressing part */
            str = meca.NP1.strike + 180.;
            while (str >= meca.NP1.strike - EPSIL) {
                i++;
                //sincosd (str, &si, &co);
                let si = Math.sin(str * D2R)
                let co = Math.cos(str * D2R)
                x[i] = x0 + si * radius_size;
                y[i] = y0 + co * radius_size;
                str -= increment;
            }
            let npoints = i + 1;
            //ps_polygon (x, y, npoints, rgb, outline);
            collect.push({x: x.slice(0,npoints), y: y.slice(0,npoints) })
            color.push(rgb)
            i = -1;
        }
        /* second nodal plane part */
        str = meca.NP2.strike;
        while (str <= meca.NP2.strike + 180. + EPSIL) {
            i++;
            let radius = proj_radius(meca.NP2.strike, meca.NP2.dip, str) * radius_size;
            //sincosd (str, &si, &co);
            let si = Math.sin(str * D2R)
            let co = Math.cos(str * D2R)
            x[i] = x0 + radius * si;
            y[i] = y0 + radius * co;
            str += increment;
        }
        if (t_axis_vertical) {
            /* inverse fault, close compressing part */
            let npoints = i+1;
            //ps_polygon(x, y, npoints, rgb, outline);
            collect.push({x: x.slice(0,npoints), y: y.slice(0,npoints) })
            color.push(rgb)
        }
        if(!t_axis_vertical) {
            /* normal fault, close second compressing part */
            str = meca.NP2.strike + 180.;
            while (str >= meca.NP2.strike - EPSIL) {
                i++;
                //sincosd (str, &si, &co);
                let si = Math.sin(str * D2R)
                let co = Math.cos(str * D2R)
                x[i] = x0 + si * radius_size;
                y[i] = y0 + co * radius_size;
                str -= increment;
            }
            let npoints = i + 1;
            collect.push({x: x.slice(0,npoints), y: y.slice(0,npoints) })
            color.push(rgb)
            //ps_polygon(x, y, npoints, rgb, outline);
        }
    }
    /* pure strike-slip */
    else if ((90. - meca.NP1.dip) < EPSIL && (90. - meca.NP2.dip) < EPSIL) {
        // pure strike slip fault
        // first compressing part
        let i = 0;
        let str = meca.NP1.strike;
        let increment = strike_slip_taxis(meca);

        while (increment > 0. ? str <= meca.NP1.strike + 90. : str >= meca.NP1.strike - 90.) {
            let si = Math.sin(str * D2R)
            let co = Math.cos(str * D2R)
            x[i] = x0 + si * radius_size;
            y[i] = y0 + co * radius_size;
            str += increment;
            i++;
        }
        x[i] = x0;
        y[i] = y0;
        let npoints = i + 1;
        collect.push({x: x.slice(0,npoints), y: y.slice(0,npoints) })
        color.push(rgb)
        //ps_polygon(x, y, npoints, rgb, outline);

        /* second compressing part */
        i = 0;
        str = meca.NP1.strike + 180.;
        while (increment > 0. ?  str <= meca.NP1.strike + 270. + EPSIL : str >= meca.NP1.strike + 90. - EPSIL) {
            let si = Math.sin(str * D2R)
            let co = Math.cos(str * D2R)
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
        let i = -1;
        let increment = 1.;
        if (meca.NP1.strike > N_axis.strike)
            meca.NP1.strike -= 360.;
        let str = meca.NP1.strike;
        log('strike p0', str,Math.abs(90. - meca.NP1.dip) < EPSIL ? meca.NP1.strike : N_axis.strike, increment )
        while (Math.abs(90. - meca.NP1.dip) < EPSIL ? str <= meca.NP1.strike + EPSIL : str <= N_axis.strike + EPSIL) {
            i++;
            let radius = proj_radius(meca.NP1.strike, meca.NP1.dip, str) * radius_size;
            let si = Math.sin(str * D2R)
            let co = Math.cos(str * D2R)
            //sincosd (str, &si, &co);
            x[i] = x0 + radius * si;
            y[i] = y0 + radius * co;
            str += increment;
        }

        /* second nodal plane from null axis */
        log('strike p1 (pre)', meca.NP2.strike, increment, fault)
        meca.NP2.strike += (1. + fault) * 90.;
        log('strike p1 (pre)', meca.NP2.strike, increment, fault)
        if (meca.NP2.strike >= 360.) meca.NP2.strike -= 360.;
        increment = fault;
        if (fault * (meca.NP2.strike - N_axis.strike) < -EPSIL) meca.NP2.strike += fault * 360.;
        str = Math.abs(90. - meca.NP2.dip) < EPSIL ? meca.NP2.strike : N_axis.strike;
        log('strike p1', str, meca.NP2.strike, increment, fault)
        while (increment > 0. ? str <= meca.NP2.strike + EPSIL : str >= meca.NP2.strike - EPSIL) {
            i++;
            let radius = proj_radius(meca.NP2.strike - (1. + fault) * 90., meca.NP2.dip, str) * radius_size;
            let si = Math.sin(str * D2R)
            let co = Math.cos(str * D2R)
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
        log('strike p2', str, meca.NP1.strike, increment)
        while (increment > 0. ? str <= meca.NP1.strike + EPSIL : str >= meca.NP1.strike - EPSIL) {
            i++;
            let si = Math.sin(str * D2R)
            let co = Math.cos(str * D2R)
            //sincosd (str, &si, &co);
            x[i] = x0 + si * radius_size;
            y[i] = y0 + co * radius_size;
            str += increment;
        }
        log('strike p3', str)

        let npoints = i + 1;
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
            let radius = proj_radius(meca.NP1.strike - 180., meca.NP1.dip, str) * radius_size;
            let si = Math.sin(str * D2R)
            let co = Math.cos(str * D2R)
            //sincosd (str, &si, &co);
            x[i] = x0 + radius * si;
            y[i] = y0 + radius * co;
            str += increment;
            //console.log('str1',str, radius, meca.NP1.dip)
        }
        log('strike p4', str)

        /* second nodal plane from null axis */
        meca.NP2.strike = zero_360(meca.NP2.strike + 180.);
        increment = -fault;
        if (fault * (N_axis.strike - meca.NP2.strike) < - EPSIL) meca.NP2.strike -= fault * 360.;
        str = Math.abs(90. - meca.NP2.dip) < EPSIL ? meca.NP2.strike : N_axis.strike;
        while (increment > 0. ? str <= meca.NP2.strike + EPSIL : str >= meca.NP2.strike - EPSIL) {
            i++;
            let radius = proj_radius(meca.NP2.strike - (1. - fault) * 90., meca.NP2.dip, str) * radius_size;
            let si = Math.sin(str * D2R)
            let co = Math.cos(str * D2R)
            //sincosd (str, &si, &co);
            x[i] = x0 + radius * si;
            y[i] = y0 + radius * co;
            str += increment;
            //console.log('str2',str)
        }
        log('strike p5', str)

        /* close the second compressing part */
        meca.NP1.strike = zero_360(meca.NP1.strike);
        meca.NP2.strike = zero_360(meca.NP2.strike);
        increment = pos_NP1_NP2 >= 0. ? -fault : fault;
        if (increment * (meca.NP1.strike - meca.NP2.strike) < - EPSIL) meca.NP1.strike += increment * 360.;
        str = meca.NP2.strike;
        while (increment > 0. ? str <= meca.NP1.strike + EPSIL : str >= meca.NP1.strike - EPSIL) {
            i++;
            let si = Math.sin(str * D2R)
            let co = Math.cos(str * D2R)
            //sincosd (str, &si, &co);
            x[i] = x0 + si * radius_size;
            y[i] = y0 + co * radius_size;
            str += increment;
            //console.log('str3',str, radius_size)
        }
        log('strike p6', str)

        npoints = i + 1;
        collect.push({x: x.slice(0,npoints), y: y.slice(0,npoints) })
        color.push(rgb)
        //ps_polygon(x, y, npoints, rgb, outline);
        {
            //console.log(meca.ax.t.strike, meca.ax.t.dip, meca.ax.t.n)
            //console.log(meca.ax.p.strike, meca.ax.p.dip, meca.ax.p.n)
            let inside1 = false
            let dip = meca.ax.t.dip
            if(close(dip, 0.0))
                dip = 1e-1
            const t = axis2xy(meca.ax.t.strike, dip)
            t.x *= radius_size
            t.y *= radius_size
            t.x += x0
            t.y += y0
            inside1 = isPointInsidePolygonRCA([t.x,t.y], collect[1].x,  collect[1].y)
            if(collect.length > 2) {
                let inside2 = isPointInsidePolygonRCA([t.x,t.y], collect[2].x,  collect[2].y)
                inside1 = inside1 || inside2
            }
            if(!inside1) {
                console.log('t inside', inside1, fault, t, meca.ax.t.strike, meca.ax.t.dip)
            }
        }
    }

    return { collect, colors: color }
}


// Formatting for number and arrays of numbers
function _f(x) {
    return x.toFixed(7).padStart(12, ' ')
}
function _fa(x) {
    return x.map(_f).join(", ")
}

function ps_tensor(x0, y0, size,  T,  N, P, c_rgb, e_rgb, outline, plot_zerotrace) {

    const D2R = Math.PI / 180.0
    const EPSIL = 1e-5
    let b = 1;
    let j = 1
    let j2 = 0
    let j3 = 0
    let lineout = 1;

    let a = zeros(3)
    let p = zeros(3)
    let v = zeros(3)
    let az = 0
    let azp = 0

    let collect = []
    let colors = []
    let azi = zeros([3,2])

    let n = 0

    let x = zeros(400)
    let y = zeros(400)
    let x2 = zeros(400)
    let y2 = zeros(400)
    let x3 = zeros(400)
    let y3 = zeros(400)

    let xp1 = zeros(800)
    let yp1 = zeros(800)
    let xp2 = zeros(400)
    let yp2 = zeros(400)

    if(Math.abs(N.val) < EPSIL && Math.abs(P.val) > EPSIL && Math.abs(T.val) > EPSIL) {
        N.val = 1e-8
    }
    log('T',T)
    log('N',N)
    log('P',P);

    a[0] = T.strike; a[1] = N.strike; a[2] = P.strike;
    p[0] = T.dip; p[1] = N.dip; p[2] = P.dip;
    v[0] = T.val; v[1] = N.val; v[2] = P.val;

    let inc = 1
    if(Math.abs(v[0]) < EPSIL && Math.abs(v[1]) < EPSIL) {
        inc = 0.1
    }
    const v0 = [...v]
    log('v: ', v[0],v[1],v[2])
    let vi = (v[0] + v[1] + v[2]) / 3.;
    for (let i=0; i<=2; i++)
        v[i] = v[i] - vi;
    log('v: ', v[0],v[1],v[2])

    if(close(len(v),0)) {
        log("EXLOSION / IMPLOSION")
        collect.push( xcircle([x0, y0], size * 0.5))
        colors.push((v0[0] > 0) ? c_rgb : e_rgb)
        return { collect, colors: colors }
    }

    /*Redundant cases for pure explosion and large isotropic with no nodes removed*/
    /*DS Dreger 2011*/

    /*Test to determine if angle exceeds plotting capability*/
    /*Used to choose the dominant eigenvalue after Frohlich for plotting purposes*/
    /*DS Dreger 2011*/
    let bigisotestv0=0
    let bigisotestv2=0
    let max_s21 = 0
    let max_s22 = 0
    let s2alphan = 0
    //console.log('possible f,iso1', -v[1]/v[0], vi/v[0],'f,iso2', -v[1]/v[2], vi/v[2])
    for (let i=0; i<360; i+=inc) {
        let fir =  i * D2R;
        let f = -v[1]/v[0];
        let iso = vi/v[0];
        s2alphan = (2. + 2. * iso) / (3. + (1. - 2. * f) * Math.cos(2. * fir));
        if(s2alphan > 1.0 )
            bigisotestv0 += 1;
        max_s21 = Math.max(max_s21, s2alphan)
        f = -v[1]/v[2];
        iso = vi/v[2];
        if(f == -1 && (i == 90  || i == 270)) {
            bigisotestv2 += 1
        } else {
            s2alphan = (2. + 2. * iso) / (3. + (1. - 2. * f) * Math.cos(2. * fir));
            //console.log(i, f, iso, s2alphan, 3 + (1-2*f) * Math.cos(2*fir), max_s22, bigisotestv2)
            if(s2alphan > 1.0 )
                bigisotestv2 += 1;
            max_s22 = Math.max(max_s22, s2alphan)
        }
    }
    /*End Test*/

    log('v: ', v[0],v[1],v[2])
    log('ps_tensor max s2', max_s21, max_s22)
    log('big iso test', bigisotestv0, bigisotestv2)

    let radius_size = size * 0.5;

    /*Determine the eigenvalue to use based on s2alphan test above*/
    /*Define the shading convention*/
    /*DS Dreger 2011*/
    let d = 0
    let m = 2
    let rgb1 = c_rgb
    let rgb2 = e_rgb
    if(bigisotestv0 == 0 ) {
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
        console.log("Error bigisotest failed", bigisotestv0, bigisotestv2)
        console.log('v: ', _fa(v), _fa(v0))
        console.log('ps_tensor max s2', max_s21, max_s22)
        console.log('big iso test', bigisotestv0, bigisotestv2)
        console.log('f1,f2    ', v[1]/v[0], v[1]/v[2])
        console.log('iso1,iso2', vi/v[0], vi/v[2])
        return undefined
    }
    log('colors',rgb1, rgb2)
    if (plot_zerotrace)
        vi = 0.;

    let f = - v[1] / v[d];
    let iso = vi / v[d];

    let spd = Math.sin(p[d] * D2R)
    let cpd = Math.cos(p[d] * D2R)
    let spb = Math.sin(p[b] * D2R)
    let cpb = Math.cos(p[b] * D2R)
    let spm = Math.sin(p[m] * D2R)
    let cpm = Math.cos(p[m] * D2R)

    let sad = Math.sin(a[d] * D2R)
    let cad = Math.cos(a[d] * D2R)
    let sab = Math.sin(a[b] * D2R)
    let cab = Math.cos(a[b] * D2R)
    let sam = Math.sin(a[m] * D2R)
    let cam = Math.cos(a[m] * D2R)

    // Cliff Frohlich, Seismological Research letters,
    // Vol 7, Number 1, January-February, 1996
    // Unless the isotropic parameter lies in the range
    // between -1 and 1 - f there will be no nodes whatsoever

    if(close(iso, 1-f)) {
        log("NO NODES")
        collect.push( xcircle([x0, y0], size * 0.5))
        colors.push((v0[0] > 0) ? c_rgb : e_rgb)
        return { collect, colors: colors }

    }
    //console.log("p", spd, spb, spm, cpd, cpb, cpm);
    //console.log("a", sad, sab, sam, cad, cab, cam);

    for (let i=0; i<360; i+=inc) {
        let fir = i * D2R;
        let s2alphan = (2. + 2. * iso) / (3. + (1. - 2. * f) * Math.cos(2. * fir));

        let alphan = Math.asin(Math.sqrt(s2alphan));
        let sfi = Math.sin(fir)
        let cfi = Math.cos(fir)
        //sincos (fir, &sfi, &cfi);
        let san = Math.sin(alphan)
        let can = Math.cos(alphan)
        //sincos (alphan, &san, &can);

        let xz = can * spd + san * sfi * spb + san * cfi * spm;
        let xn = can * cpd * cad + san * sfi * cpb * cab + san * cfi * cpm * cam;
        let xe = can * cpd * sad + san * sfi * cpb * sab + san * cfi * cpm * sam;

        let takeoff = 0
        az = 0

        if (Math.abs(xn) < EPSIL && Math.abs(xe) < EPSIL) {
            takeoff = 0.;
            az = 0.;
        } else {
            az = Math.atan2(xe, xn);
            if (az < 0.)
                az += Math.PI * 2.;
            takeoff = Math.acos(xz / Math.sqrt(xz * xz + xn * xn + xe * xe));
        }
        if (takeoff >= (Math.PI / 2.0)) {
            takeoff = Math.PI - takeoff;
            az += Math.PI;
            if (az > Math.PI * 2.)
                az -= Math.PI * 2.;
        }

        //console.log(i, n, az / D2R, takeoff)
        let r = Math.sqrt(2.0) * Math.sin(takeoff / 2.);
        let si = Math.sin(az)
        let co = Math.cos(az)
        //console.log(i, n, _f(xe), _f(xn), _f(xz), _f(az / D2R), _f(takeoff / D2R), _f(r))
        //sincos (az, &si, &co);
        if (i == 0) {
            azi[i][0] = az;
            x[i] = x0 + radius_size * r * si;
            y[i] = y0 + radius_size * r * co;
            azp = az;
        } else {
            if (Math.abs(Math.abs(az - azp) - Math.PI) < D2R * 10. && takeoff > 80 * D2R) {
                //console.log('*******n',i, n,n+1, az/D2R, azp/D2R, Math.abs(Math.abs(az - azp) - Math.PI), 10.0 * D2R)
                azi[n][1] = azp;
                n += 1
                if(n >= 3) {
                    if(Math.abs(N.val) < 1e-9) {
                        console.error('More than three sections', n, N.val)
                        N.val = -1e-8
                        return ps_tensor(x0, y0, size, T, N, P, c_rgb, e_rgb, outline, plot_zerotrace)
                    }
                    console.error('More than three sections', n, N.val)
                    return undefined;
                    //return ps_tensor_n(x0, y0, size, T, N, P, c_rgb, e_rgb, outline, plot_zerotrace)
                } else {
                    azi[n][0] = az;
                }
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
    if(n < 3)
        azi[n][1] = az;

    //console.log('n',n)
    //console.log('azi',azi.map(v => [v[0]/D2R, v[1]/D2R]))
    //console.log('color', rgb2)
    collect.push( xcircle([x0, y0], radius_size ))
    colors.push(rgb2)
    //ps_circle(x0, y0, radius_size*2., rgb2, lineout);

    for(let i = 0; i < n+1; i++) {
        log(i, azi[i][0]/D2R, azi[i][1] / D2R)
    }
    if(n > 3) {
        return undefined;
    }
    if(n == 0) {
        let npoints = 360;
        //ps_polygon(xp1, yp1, npoints, rgb1, outline);
        collect.push({ x: x.slice(0, npoints), y: y.slice(0, npoints) })
        colors.push(rgb1)
    }
    if(n == 1) {
        //case 1 :
        let i = 0
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
        let npoints = i;
        //ps_polygon(xp1, yp1, npoints, rgb1, outline);
        collect.push({ x: xp1.slice(0, npoints), y: yp1.slice(0, npoints) })
        colors.push(rgb1)

        x = []
        y = []
        i = 0
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
        for (let ii=0; ii<j; ii++) {
            xp1[i] = x[ii];
            yp1[i] = y[ii];
            i += 1
        }
        if (azi[2][0] - azi[0][1] > Math.PI)
            azi[2][0] -= Math.PI * 2.;
        else if (azi[0][1] - azi[2][0] > Math.PI)
            azi[2][0] += Math.PI * 2.;
        if (azi[2][0] < azi[0][1]) {
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

        let npoints = i;
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

export function plane_normal(strike, dip) {
    const D2R = Math.PI / 180.0

    // Right Hand Coordinate system
    // x - Right
    // y - Forward (90 degrees counter clockwise from X)
    // z - Up

    // Convert from N-clockwise (Geographic) to x-y-z coordinates
    // With a rotation for trailing normal if pointing downward
    strike = 90 - strike - 90

    // Planes should have a downward pointing normal
    // Adding 90 for co-latitude and -90 to get a downward pointing normal
    dip = dip + 90 + 90

    const x =  Math.sin(dip * D2R) * Math.cos(strike * D2R)
    const y =  Math.sin(dip * D2R) * Math.sin(strike * D2R)
    const z =  Math.cos(dip * D2R)

    return norm([x,y,z])
}
export function normal_plane(n0) {
    // Computes strike and dip from Plane Normal
    const EPSIL = 1e-5
    const R2D = 180.0 / Math.PI
    const D2R = Math.PI / 180.0
    let n = norm(n0)
    n = downward_normal(n)
    const r = len(n)
    if(close(n[2],1) || close(n[2],-1)) {
        return { str: 0, dip: 0  }
    }
    const theta = R2D * Math.acos(n[2] / r)
    const phi = R2D * Math.atan2(n[1],n[0])
    log('normal_plane (phi/strike)',phi, n[1],n[0])

    const str = zero_360(90 - phi + 90)
    const dip = 90 + 90 - theta
    return { str, dip }
}

function normal_to_strike_vec(n) {
    const EPSIL = 1e-5
    const vec = norm(n)
    // Horizontal plane "defined" with a strike of 0 deg.
    if(Math.abs(vec[2] - 1.0) < EPSIL)
        return [0, 1, 0]
    if(Math.abs(vec[2] - - 1.0) < EPSIL)
        return [0, 1, 0]

    // Positive rotation from normal, +90 degrees clockwise, leading
    return norm([ vec[1], -vec[0], 0])
}

// A * scalar
function fM(a, fac) {
    let m = [[0,0,0],[0,0,0],[0,0,0]]
    for(let i = 0; i < 3; i++) {
        for(let j = 0; j < 3; j++) {
            m[i][j] = a[i][j] * fac;
        }
    }
    return m
}

// A + B
function addM(a,b) {
    let m = [[0,0,0],[0,0,0],[0,0,0]]
    for(let i = 0; i < 3; i++) {
        for(let j = 0; j < 3; j++) {
            m[i][j] = a[i][j] + b[i][j];
        }
    }
    return m
}
// A * vectpr
function mulM(a,v) {
    let out = [0,0,0]
    for(let i = 0; i < 3; i++) {
        for(let j= 0; j < 3; j++) {
            out[i] += a[i][j] * v[j]
        }
    }
    return out
}

export function normal_second_plane(n, rake, strike) {
    // Computes Normal to the second plane from Normal to first plane and rake
    const EPSIL = 1e-5
    const D2R = Math.PI / 180.0
    let str_vec = normal_to_strike_vec(n)
    if(close(n[2],1) || close(n[2],-1)) { // Horizontal Plane, strike defined by strike value
        str_vec = [ Math.cos((90-strike) * D2R), Math.sin((90-strike) * D2R) , 0]
    }
    const rake_r = -rake * D2R
    // https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
    //   R = ( cos θ ) I + ( sin θ ) [ u ] × + ( 1 − cos⁡ θ ) ( u ⊗ u )
    const ux = [[   0,  -n[2],  n[1]],
                [ n[2],    0,  -n[0]],
                [-n[1],  n[0],     0]]
    const uxu = [[n[0]*n[0], n[1]*n[0], n[2]*n[0]],
                 [n[0]*n[1], n[1]*n[1], n[2]*n[1]],
                 [n[0]*n[2], n[1]*n[2], n[2]*n[2]]]

    const I = [[1,0,0],[0,1,0],[0,0,1]]
    let R = addM(fM(I, Math.cos(rake_r)), fM(ux, Math.sin(rake_r)))
    R = addM(R, fM(uxu, 1.0 - Math.cos(rake_r)))
    return mulM(R, str_vec)
}

function normals_to_axes(n1, n2) {
    const ax1 = [n1[0]+n2[0], n1[1]+n2[1], n1[2]+n2[2]]
    const ax2 = [n1[0]-n2[0], n1[1]-n2[1], n1[2]-n2[2]]
    return { ax1, ax2 }
}

function axes_to_planes_dc(n1, n2) {
    // Convert P and T axes to Planes (Strike, dip, rake) + Plane normals
    //console.log('n1',_fa(n1.v))
    //console.log('n2',_fa(n2.v))
    const ax = normals_to_axes(n1.v, n2.v)
    if(ax.ax1[2] > 1e-5)
        ax.ax1 = ax.ax1.map(v => -v)
    if(ax.ax2[2] > 1e-5)
        ax.ax2 = ax.ax2.map(v => -v)

    const pole1 = new Pole(ax.ax1)
    const pole2 = new Pole(ax.ax2)
    const plane1 = pole1.plane()
    const plane2 = pole2.plane()
    //console.log('plane1', plane1)
    //console.log('plane2', plane2)

    //console.log("calling rake3 from axes_to_planes_dc")
    const rakes = rake3(pole1.n, pole2.n, n1.v, n2.v, plane1.strike, plane2.strike)

    const np1 = pole1.plane().to_nodal_plane(rakes.rake1)
    const np2 = pole2.plane().to_nodal_plane(rakes.rake2)
    //console.log('np1', np1)
    //console.log('np2', np2)
    return {NP1: np1, NP2: np2}
}

function downward_normal(n) {
    const EPSIL = 1e-5
    const fac = (n[2] > EPSIL) ? -1 : 1
    return n.map(v => fac * v)
}

export function auxillary_plane(strike, dip, rake) {
    const np = new NodalPlane(strike, dip, rake)
    const aux = np.aux_plane()
    aux.str = aux.strike
    aux.rake0 = rake
    return aux
}

export function dc2axe (meca) { //, struct SEIS_AXIS *T, struct SEIS_AXIS *N, struct SEIS_AXIS *P) {
  /*
    From FORTRAN routines of Anne Deschamps :
    compute azimuth and plungement of P-T axis
    from nodal plane strikes, dips and rakes.
  */
    const D2R = Math.PI / 180.0
    const R2D = 180.0 / Math.PI
    const M_SQRT2 = Math.sqrt(2.0)
    const EPSIL = 1e-5
    //double cd1, sd1, cd2, sd2, cp1, sp1, cp2, sp2;
    //double amz, amx, amy, dx, px, dy, py;

    let cd1 = Math.cos (meca.NP1.dip * D2R) * M_SQRT2
    let sd1 = Math.sin (meca.NP1.dip * D2R) * M_SQRT2
    let cd2 = Math.cos (meca.NP2.dip * D2R) * M_SQRT2
    let sd2 = Math.sin (meca.NP2.dip * D2R) * M_SQRT2
    let cp1 = - Math.cos (meca.NP1.strike * D2R) * sd1; // r sin theta cos phi = y
    let sp1 = Math.sin (meca.NP1.strike * D2R) * sd1;   // r sin theta sin phi = x
    let cp2 = - Math.cos (meca.NP2.strike * D2R) * sd2;
    let sp2 = Math.sin (meca.NP2.strike * D2R) * sd2;

    let amz = - (cd1 + cd2);
    let amx = - (sp1 + sp2);
    let amy = cp1 + cp2;
    let N1 = [-amx, amy, amz]
    let dx = Math.atan2 (Math.sqrt(amx*amx + amy*amy), amz) * R2D - 90.0;
    let px = Math.atan2 (amy, -amx) * R2D;
    if (px < 0.0)
        px += 360.0;
    if (dx < EPSIL) {
        if (px > 90.0 && px < 180.0)
            px += 180.0;
        if (px >= 180.0 && px < 270.0)
            px -= 180.0;
    }

    amz = cd1 - cd2;
    amx = sp1 - sp2;
    amy = - cp1 + cp2;
    let N2 = [amx, amy, amz]
    let dy = Math.atan2 (Math.sqrt(amx*amx + amy*amy), -Math.abs(amz)) * R2D - 90.0;
    let py = Math.atan2 (amy, -amx) * R2D;
    if (amz > 0.0)
        py -= 180.0;
    if (py < 0.0)
        py += 360.0;
    if (dy < EPSIL) {
        if (py > 90.0 && py < 180.0)
            py += 180.0;
        if (py >= 180.0 && py < 270.0)
            py -= 180.0;
    }
    let P = {strike: 0, dip: 0}
    let T = {strike: 0, dip: 0}
    let N = {str: 0, dip: 0}
    if (meca.NP1.rake > 0.0) {
        P.dip = dy; P.strike = py;
        T.dip = dx; T.strike = px;
    } else {
        P.dip = dx; P.strike = px;
        T.dip = dy; T.strike = py;
    }

    N.str = null_axis_strike (T.strike, T.dip, P.strike, P.dip);
    N.dip = null_axis_dip (T.strike, T.dip, P.strike, P.dip);

    T.n = norm(N1)
    P.n = norm(N2)

    return {n: N, p: P, t: T}
}


// Check if a point is within a polygon (pts)
//  https://en.wikipedia.org/wiki/Point_in_polygon#Ray_casting_algorithm
function isPointInsidePolygonRCA(point, xv, yv) {
    let n = xv.length;
    let xp = point[0];
    let yp = point[1];
    //let xv = pts.map(p => p[0]);
    //let yv = pts.map(p => p[1]);

  if (Math.abs(xv[0] - xv[n - 1]) < 1e-7 && Math.abs(yv[0] - yv[n - 1]) < 1e-7) {
    n -= 1;
  }
  let x2 = xv[n - 1]
  let y2 = yv[n - 1]
  let nleft = 0

  let x1 = x2;
  let y1 = y2;

  // Loop over line segments (assuming the polygon is closed)
  for (let i = 0; i < n; i++) {
    x1 = x2
    y1 = y2
    x2 = xv[i]
    y2 = yv[i]
    if (y1 >= yp && y2 >= yp) {
      continue;
    }
    if (y1 < yp && y2 < yp) {
      continue;
    }
    if (y1 == y2) {
      if (x1 >= xp && x2 >= xp) {
        continue;
      }
      if (x1 < xp && x2 < xp) {
        continue;
      }
      nleft += 1;
    } else {
      let xi = x1 + (yp - y1) * (x2 - x1) / (y2 - y1);
      if (xi == xp) {
        nleft = 1;
        break;
      }
      if (xi > xp) {
        nleft += 1;
      }
    }
  }
  let xin = nleft % 2;
  return xin == 1;
}

export function axis2xy(strike, dip) {
    const D2R = Math.PI / 180.0
    let radius = Math.sqrt(1. - Math.sin(dip * D2R));
    let x = radius * Math.sin(strike * D2R)
    let y = radius * Math.cos(strike * D2R)
    return {x, y}
}
