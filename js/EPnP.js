import {Matrix, QR} from 'ml-matrix';
import {SVD} from 'svd-js'

function dot(vec0, vec1) {
  return vec0[0] * vec1[0] + vec0[1] * vec1[1] + vec0[2] * vec1[2];
}

function distSquare(p0, p1) {
  return (p0[0] - p1[0]) * (p0[0] - p1[0]) +
    (p0[1] - p1[1]) * (p0[1] - p1[1]) +
    (p0[2] - p1[2]) * (p0[2] - p1[2]);
}

function sortUResult(u, q) {
  let Q = new Array(q.length);
  for (let k = 0; k < Q.length; ++k) {
    Q[k] = [q[k], k];
  }
  Q.sort( (x, y) => {return y[0] - x[0]});
  let sortedQ = new Array(q.length);
  for (let k = 0; k < sortedQ.length; ++k) {
    sortedQ[k] = Q[k][0];
  }
  let sortedU = new Array(u.length);
  for (let i = 0; i < sortedU.length; ++i) {
    sortedU[i] = new Array(u[0].length);
  }
  for (let k = 0; k < q.length; ++k) {
    let oriK = Q[k][1];
    for (let i = 0; i < sortedU.length; ++i)
      sortedU[i][k] = u[i][oriK];
  }
  return {sortedU, sortedQ};
}

function svdInverse(A) {
  let {u, v, q} = SVD(A);
  for (let i = 0; i < v.length; ++i) {
    for (let j = 0; j < v[0].length; ++j) {
      v[i][j] /= q[j];
    }
  }
  let AInv = new Array(A[0].length);
  for (let i = 0; i < AInv.length; ++i) {
    AInv[i] = new Array(A.length);
    for (let j = 0; j < A.length; ++j) {
      AInv[i][j] = 0;
      for (let k = 0; k < u[0].length; ++k) {
        AInv[i][j] += v[i][k] * u[j][k];
      }
    }
  }
  return AInv;
}

function svdSolve(A, b) {
  const AInv = svdInverse(A);
  let x = new Array(A[0].length);
  for (let i = 0; i < x.length; ++i) {
    x[i] = 0;
    for (let j = 0; j < b.length; ++j) {
      x[i] += AInv[i][j] * b[j];
    }
  }
  return x;
}

class EPnPSolver {
  constructor(ptsNum, intrinsics) {
    this.pointsNumber = ptsNum;
    this.pointsNumberInverse = 1 / this.pointsNumber;
    this.fx = intrinsics[0];
    this.fy = intrinsics[1];
    this.cx = intrinsics[2];
    this.cy = intrinsics[3];
    this._initMatrices();
  }

  _initMatrices() {
    this.alphas = new Array(this.pointsNumber);
    for (let i = 0; i < this.alphas.length; ++i) {
      this.alphas[i] = [0, 0, 0, 0];
    }
    this.M = new Array(2 * this.pointsNumber);
    for (let i = 0; i < this.M.length; ++i) {
      this.M[i] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    }
    this.MtM = new Array(12);
    for (let i = 0; i < 12; ++i) {
      this.MtM[i] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    }
    this.l6x10 = new Array(6);
    for (let i = 0; i < 6; ++i) {
      this.l6x10[i] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    }
    this.l6x4 = new Array(6);
    for (let i = 0; i < 6; ++i) {
      this.l6x4[i] = [0, 0, 0, 0];
    }
    this.l6x3 = new Array(6);
    for (let i = 0; i < 6; ++i) {
      this.l6x3[i] = [0, 0, 0];
    }
    this.l6x5 = new Array(6);
    for (let i = 0; i < 6; ++i) {
      this.l6x5[i] = [0, 0, 0, 0, 0];
    }

    this.dv = [
      [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]],
      [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]],
      [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]],
      [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]]
    ];
    this.rho = [0, 0, 0, 0, 0, 0];

    this.pcs = new Array(this.pointsNumber);
    for (let i = 0; i < this.pointsNumber; ++i) {
      this.pcs[i] = new Array(3);
    }

    // coordinates in wolrd control points;
    this.pcws = new Array(this.pointsNumber);
    for (let i = 0; i < this.pointsNumber; ++i) {
      this.pcws[i] = new Array(3);
    }
  }

  _chooseControlPoints(pws) {
    let cws = [ [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (let i = 0; i < this.pointsNumber; ++i) {
      cws[0][0] += pws[i][0];
      cws[0][1] += pws[i][1];
      cws[0][2] += pws[i][2];
    }
    cws[0][0] *= this.pointsNumberInverse;
    cws[0][1] *= this.pointsNumberInverse;
    cws[0][2] *= this.pointsNumberInverse;
    
    for (let i = 0; i < this.pointsNumber; ++i) {
      this.pcws[i][0] = pws[i][0] - cws[0][0];
      this.pcws[i][1] = pws[i][1] - cws[0][1];
      this.pcws[i][2] = pws[i][2] - cws[0][2];
    }
    let pcwstpcws = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    for (let i = 0; i < 3; ++i) {
      for (let j = 0; j < 3; ++j) {
        for (let k = 0; k < this.pointsNumber; ++k) {
          pcwstpcws[i][j] += this.pcws[k][i] * this.pcws[k][j];
        }
      }
    }
    const {u, v, q} = SVD(pcwstpcws, true, false);
    const {sortedU, sortedQ} = sortUResult(u, q);
    for (let i = 1; i < 4; ++i) {
      let k = Math.sqrt(sortedQ[i-1] / this.pointsNumber);
      for (let j = 0; j < 3; ++j) {
        cws[i][j] = cws[0][j] + k * sortedU[j][i - 1];
      }    
    }
    return cws;
  }

  _computeBarycentricCoordinates(cws, pws) {
    let cc = [
      [cws[1][0] - cws[0][0], cws[2][0] - cws[0][0], cws[3][0] - cws[0][0]],
      [cws[1][1] - cws[0][1], cws[2][1] - cws[0][1], cws[3][1] - cws[0][1]],
      [cws[1][2] - cws[0][2], cws[2][2] - cws[0][2], cws[3][2] - cws[0][2]]
    ];
    let cc_inv = svdInverse(cc)
    for (let i = 0; i < this.pointsNumber; ++i) {
      this.alphas[i][1] =
        cc_inv[0][0] * this.pcws[i][0] + cc_inv[0][1] * this.pcws[i][1] + cc_inv[0][2] * this.pcws[i][2];
      this.alphas[i][2] =
        cc_inv[1][0] * this.pcws[i][0] + cc_inv[1][1] * this.pcws[i][1] + cc_inv[1][2] * this.pcws[i][2];
      this.alphas[i][3] =
        cc_inv[2][0] * this.pcws[i][0] + cc_inv[2][1] * this.pcws[i][1] + cc_inv[2][2] * this.pcws[i][2];
      this.alphas[i][0] = 1 - this.alphas[i][1] - this.alphas[i][2] - this.alphas[i][3];
    }
  }

  _fillM(pis) {
    for (let i = (this.pointsNumber - 1) >>> 0; i >= 0; --i) {
      this.M[2 * i][0] = this.alphas[i][0] * this.fx;
      this.M[2 * i][1] = 0;
      this.M[2 * i][2] = this.alphas[i][0] * (this.cx - pis[i][0]);
      this.M[2 * i][3] = this.alphas[i][1] * this.fx;
      this.M[2 * i][4] = 0;
      this.M[2 * i][5] = this.alphas[i][1] * (this.cx - pis[i][0]);
      this.M[2 * i][6] = this.alphas[i][2] * this.fx;
      this.M[2 * i][7] = 0;
      this.M[2 * i][8] = this.alphas[i][2] * (this.cx - pis[i][0]);
      this.M[2 * i][9] = this.alphas[i][3] * this.fx;
      this.M[2 * i][10] = 0;
      this.M[2 * i][11] = this.alphas[i][3] * (this.cx - pis[i][0]);
  
      this.M[2 * i + 1][0] = 0;
      this.M[2 * i + 1][1] = this.alphas[i][0] * this.fy;
      this.M[2 * i + 1][2] = this.alphas[i][0] * (this.cy - pis[i][1]);
      this.M[2 * i + 1][3] = 0;
      this.M[2 * i + 1][4] = this.alphas[i][1] * this.fy;
      this.M[2 * i + 1][5] = this.alphas[i][1] * (this.cy - pis[i][1]);
      this.M[2 * i + 1][6] = 0;
      this.M[2 * i + 1][7] = this.alphas[i][2] * this.fy;
      this.M[2 * i + 1][8] = this.alphas[i][2] * (this.cy - pis[i][1]);
      this.M[2 * i + 1][9] = 0;
      this.M[2 * i + 1][10] = this.alphas[i][3] * this.fy;
      this.M[2 * i + 1][11] = this.alphas[i][3] * (this.cy - pis[i][1]);
    }
  }

  _computeMtM() {
    for (let i = 0; i < 12; ++i) {
      for (let j = 0; j < 12; ++j) {
        this.MtM[i][j] = 0;
        for (let k = 0; k < 2 * this.pointsNumber; ++k)
          this.MtM[i][j] += this.M[k][i] * this.M[k][j];
      }
    }
  }

  _computeL6x10(u) {
    for (let i = 0; i < 4; ++i) {
      let a = 0;
      let b = 1;
      for (let j = 0; j < 6; ++j) {
        this.dv[i][j][0] = u[3 * a    ][11 - i] - u[3 * b][11 - i];
        this.dv[i][j][1] = u[3 * a + 1][11 - i] - u[3 * b + 1][11 - i];
        this.dv[i][j][2] = u[3 * a + 2][11 - i] - u[3 * b + 2][11 - i];
        b++;
        if (b > 3) {
          a++;
          b = a + 1;
        }
      }
    }
    for (let i = 0; i < 6; ++i) {
      this.l6x10[i][0] =     dot(this.dv[0][i], this.dv[0][i]);
      this.l6x10[i][1] = 2 * dot(this.dv[0][i], this.dv[1][i]);
      this.l6x10[i][2] =     dot(this.dv[1][i], this.dv[1][i]);
      this.l6x10[i][3] = 2 * dot(this.dv[0][i], this.dv[2][i]);
      this.l6x10[i][4] = 2 * dot(this.dv[1][i], this.dv[2][i]);
      this.l6x10[i][5] =     dot(this.dv[2][i], this.dv[2][i]);
      this.l6x10[i][6] = 2 * dot(this.dv[0][i], this.dv[3][i]);
      this.l6x10[i][7] = 2 * dot(this.dv[1][i], this.dv[3][i]);
      this.l6x10[i][8] = 2 * dot(this.dv[2][i], this.dv[3][i]);
      this.l6x10[i][9] =     dot(this.dv[3][i], this.dv[3][i]);
    }
  }

  _computeRho(cws) {
    this.rho[0] = distSquare(cws[0], cws[1]);
    this.rho[1] = distSquare(cws[0], cws[2]);
    this.rho[2] = distSquare(cws[0], cws[3]);
    this.rho[3] = distSquare(cws[1], cws[2]);
    this.rho[4] = distSquare(cws[1], cws[3]);
    this.rho[5] = distSquare(cws[2], cws[3]);
  }

  _findBetasApprox1(betas) {
    for (let i = 0; i < 6; ++i) {
      this.l6x4[i][0] = this.l6x10[i][0];
      this.l6x4[i][1] = this.l6x10[i][1];
      this.l6x4[i][2] = this.l6x10[i][3];
      this.l6x4[i][3] = this.l6x10[i][6];
    }
    let b4 = svdSolve(this.l6x4, this.rho);
    if (b4[0] < 0) {
      betas[0] = Math.sqrt(-b4[0]);
      betas[1] = -b4[1] / betas[0];
      betas[2] = -b4[2] / betas[0];
      betas[3] = -b4[3] / betas[0];
    } else {
      betas[0] = Math.sqrt(b4[0]);
      betas[1] = b4[1] / betas[0];
      betas[2] = b4[2] / betas[0];
      betas[3] = b4[3] / betas[0];
    }
  }

  _findBetasApprox2(betas) {
    for (let i = 0; i < 6; ++i) {
      this.l6x3[i][0] = this.l6x10[i][0];
      this.l6x3[i][1] = this.l6x10[i][1];
      this.l6x3[i][2] = this.l6x10[i][2];
    }
    let b3 = svdSolve(this.l6x3, this.rho);
    if (b3[0] < 0) {
      betas[0] = Math.sqrt(-b3[0]);
      betas[1] = (b3[2] < 0) ? Math.sqrt(-b3[2]) : 0;
    } else {
      betas[0] = Math.sqrt(b3[0]);
      betas[1] = (b3[2] > 0) ? Math.sqrt(b3[2]) : 0;
    }
    if (b3[1] < 0) betas[0] = -betas[0];
    betas[2] = 0;
    betas[3] = 0;
  }

  _findBetasApprox3(betas) {
    for (let i = 0; i < 6; ++i) {
      this.l6x5[i][0] = this.l6x10[i][0];
      this.l6x5[i][1] = this.l6x10[i][1];
      this.l6x5[i][2] = this.l6x10[i][2];
      this.l6x5[i][3] = this.l6x10[i][3];
      this.l6x5[i][4] = this.l6x10[i][4];
    }
    let b5 = svdSolve(this.l6x5, this.rho);
    if (b5[0] < 0) {
      betas[0] = Math.sqrt(-b5[0]);
      betas[1] = (b5[2] < 0) ? Math.sqrt(-b5[2]) : 0.0;
    } else {
      betas[0] = Math.sqrt(b5[0]);
      betas[1] = (b5[2] > 0) ? Math.sqrt(b5[2]) : 0.0;
    }
    if (b5[1] < 0) betas[0] = -betas[0];
    betas[2] = b5[3] / betas[0];
    betas[3] = 0.0;
  }

  _gaussNewton(betas) {
    const interations_number = 5;
    let A = new Matrix(6, 4);
    let b = Matrix.columnVector([0, 0, 0, 0, 0, 0]);
    for (let k = 0; k < interations_number; ++k) {
      for (let i = 0; i < 6; ++i) {
        A.set(i, 0, 2 * this.l6x10[i][0] * betas[0] +     this.l6x10[i][1] * betas[1] +     this.l6x10[i][3] * betas[2] +     this.l6x10[i][6] * betas[3]);
        A.set(i, 1,     this.l6x10[i][1] * betas[0] + 2 * this.l6x10[i][2] * betas[1] +     this.l6x10[i][4] * betas[2] +     this.l6x10[i][7] * betas[3]);
        A.set(i, 2,     this.l6x10[i][3] * betas[0] +     this.l6x10[i][4] * betas[1] + 2 * this.l6x10[i][5] * betas[2] +     this.l6x10[i][8] * betas[3]);
        A.set(i, 3,     this.l6x10[i][6] * betas[0] +     this.l6x10[i][7] * betas[1] +     this.l6x10[i][8] * betas[2] + 2 * this.l6x10[i][9] * betas[3]);
        
        b.set(i, 0, this.rho[i] - (
          this.l6x10[i][0] * betas[0] * betas[0] +
          this.l6x10[i][1] * betas[0] * betas[1] +
          this.l6x10[i][2] * betas[1] * betas[1] +
          this.l6x10[i][3] * betas[0] * betas[2] +
          this.l6x10[i][4] * betas[1] * betas[2] +
          this.l6x10[i][5] * betas[2] * betas[2] +
          this.l6x10[i][6] * betas[0] * betas[3] +
          this.l6x10[i][7] * betas[1] * betas[3] +
          this.l6x10[i][8] * betas[2] * betas[3] +
          this.l6x10[i][9] * betas[3] * betas[3]
        ));
      }
      let qr = new QR(A);
      let x = qr.solve(b).to1DArray();
      betas[0] += x[0];
      betas[1] += x[1];
      betas[2] += x[2];
      betas[3] += x[3];
    }
  }
  
  _computeRT(pws, u, betas) {
    // compute control points in camera frame
    let ccs = [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (let i = 0; i < 4; ++i) {
      for (let j = 0; j < 4; ++j) {
        for (let k = 0; k < 3; ++k) {
          ccs[j][k] += betas[i] * u[3 * j + k][11 - i];
        }
      }
    }
    // compute points in camera frame
    for (let i = 0; i < this.pointsNumber; ++i) {
      for (let j = 0; j < 3; ++j) {
        this.pcs[i][j] = this.alphas[i][0] * ccs[0][j] +
          this.alphas[i][1] * ccs[1][j] +
          this.alphas[i][2] * ccs[2][j] +
          this.alphas[i][3] * ccs[3][j];
      }
    }
    // if z less than zero
    if (this.pcs[0][2] < 0.0) {
      for(let i = 0; i < 4; i++)
        for(let j = 0; j < 3; j++)
          ccs[i][j] = -ccs[i][j];
  
      for(let i = 0; i < this.pointsNumber; i++) {
        this.pcs[i][0] *= -1;
        this.pcs[i][1] *= -1
        this.pcs[i][2] *= -1;
      }
    }

    // estimate r and t
    let pc0 = [0, 0, 0];
    let pw0 = [0, 0, 0]
    for (let i = 0; i < this.pointsNumber; ++i) {
      for (let j = 0; j < 3; ++j) {
        pc0[j] += this.pcs[i][j];
        pw0[j] += pws[i][j];
      }
    }
    for(let j = 0; j < 3; ++j) {
      pc0[j] *= this.pointsNumberInverse;
      pw0[j] *= this.pointsNumberInverse;
    }

    let abt = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    for(let i = 0; i < this.pointsNumber; i++) {
      for(let j = 0; j < 3; j++) {
        abt[j][0] += (this.pcs[i][j] - pc0[j]) * (pws[i][0] - pw0[0]);
        abt[j][1] += (this.pcs[i][j] - pc0[j]) * (pws[i][1] - pw0[1]);
        abt[j][2] += (this.pcs[i][j] - pc0[j]) * (pws[i][2] - pw0[2]);
      }
    }
    const {u: abtU, v: abtV, q} = SVD(abt);
    let R = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (let i = 0; i < 3; ++i) {
      for (let j = 0; j < 3; ++j) {
        R[i][j] = dot(abtU[i], abtV[j]);
      }
    }
    const det =
      R[0][0] * R[1][1] * R[2][2] + R[0][1] * R[1][2] * R[2][0] + R[0][2] * R[1][0] * R[2][1] -
      R[0][2] * R[1][1] * R[2][0] - R[0][1] * R[1][0] * R[2][2] - R[0][0] * R[1][2] * R[2][1];
    if (det < 0) {
      R[2][0] = -R[2][0];
      R[2][1] = -R[2][1];
      R[2][2] = -R[2][2];
    }
    let T = [
      pc0[0] - dot(R[0], pw0),
      pc0[1] - dot(R[1], pw0),
      pc0[2] - dot(R[2], pw0),
    ];
    return {R, T};
  }

  _computeReprojectionError(pws, pis, R, T) {
    let sum = 0;
    for (let i = 0; i < this.pointsNumber; ++i) {
      let x = dot(R[0], pws[i]) + T[0];
      let y = dot(R[1], pws[i]) + T[1];
      let zInv = 1 / (dot(R[2], pws[i]) + T[2]);
      let u = this.cx + this.fx * x * zInv;
      let v = this.cy + this.fy * y * zInv;

      sum += Math.sqrt((pis[i][0] - u) * (pis[i][0] - u) + (pis[i][1] - v) * (pis[i][1] - v));
    }
    return sum * this.pointsNumberInverse;
  }

  solvePnP(pis, pws) {
    let cws = this._chooseControlPoints(pws);
    this._computeBarycentricCoordinates(cws, pws);
    this._fillM(pis);
    this._computeMtM();
    const {u, v, q} = SVD(this.MtM, true, false);
    const {sortedU, sortedQ} = sortUResult(u, q);
    this._computeL6x10(sortedU);
    this._computeRho(cws);
    
    let betas = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]];
    this._findBetasApprox1(betas[1]);
    this._gaussNewton(betas[1]);
    let {R: R1, T: T1} = this._computeRT(pws, sortedU, betas[1]);
    let err1 = this._computeReprojectionError(pws, pis, R1, T1);

    this._findBetasApprox2(betas[2]);
    this._gaussNewton(betas[2]);
    let {R: R2, T: T2} = this._computeRT(pws, sortedU, betas[2]);
    let err2 = this._computeReprojectionError(pws, pis, R2, T2);

    this._findBetasApprox3(betas[3]);
    this._gaussNewton(betas[3]);
    let {R: R3, T: T3} = this._computeRT(pws, sortedU, betas[3]);
    let err3 = this._computeReprojectionError(pws, pis, R3, T3);
    let R = R1;
    let T = T1;
    if (err2 < err1) {
      if (err2 < err3) {
        R = R2;
        T = T2;
      } else {
        R = R3;
        T = T3;
      }
    } else {
      if (err3 < err1) {
        R = R3;
        T = T3;
      }
    }
    return {R, T};
  }
}

export { EPnPSolver as default };