import matrix from 'matrix-js';
import {SVD} from "svd-js";
import {Vector3} from 'three';


function procrustes(worldLandmarks) {
    // copy in this method is so stupid!!!!!!!!

    const inputLandmarksArray = [];
    const inputLandmarksArrayCopy = []
    worldLandmarks.forEach((value) => {
        inputLandmarksArray.push(
            [value.x, value.y, value.z]
        );
        inputLandmarksArrayCopy.push(
            [value.x, value.y, value.z]
        );
    });
    // const inputLandmarksArrayCopy = Object.values(inputLandmarksArray);

    const vertices = [
        [15.163, 5.0267, 143.78],
        [-15.163, 5.0267, 143.78],
        [8.2078, -0.47654, 97.24],
        [-8.2078, -0.47653, 97.24]
    ];
    const verticesCopy = [
        [15.163, 5.0267, 143.78],
        [-15.163, 5.0267, 143.78],
        [8.2078, -0.47654, 97.24],
        [-8.2078, -0.47653, 97.24]
    ];
    const verticesCopy2 = [
        [15.163, 5.0267, 143.78],
        [-15.163, 5.0267, 143.78],
        [8.2078, -0.47654, 97.24],
        [-8.2078, -0.47653, 97.24]
    ];

    console.log("input", inputLandmarksArrayCopy, verticesCopy2, vertices);

    // procrustes test input
    // inputLandmarksArray = [
    //     [0.28053278, 0.88696484, 0],
    //     [-0.28053278, 0.88696484, 0],
    //     [0.16824273, 0, 0],
    //     [-0.16824273, 0, 0]
    // ]
    //
    // vertices = [
    //     [0.09714, 0.021801, 0.78685],
    //     [-0.09714, 0.021801, 0.78685],
    //     [0.050814, -0.006826, 0.52545],
    //     [-0.050814, -0.006826, 0.52545]
    // ]

    // procrustes test output should be closed to:
    // rotation: [[ 1.00000000e+00 -2.67088823e-18 -5.25079952e-18]
    //          [-4.92883095e-18  1.08863283e-01 -9.94056732e-01]
    //          [ 3.22663370e-18  9.94056732e-01  1.08863283e-01]]
    // translation : [-6.81703893e-18 -1.69663065e+00 -2.09689233e-01]

    // solve center
    const center1 = new Vector3(0, 0, 0);
    const center2 = new Vector3(0, 0, 0);

    inputLandmarksArray.forEach((p) => {
        center1.x += p[0];
        center1.y += p[1];
        center1.z += p[2];
    });
    verticesCopy.forEach((p) => {
        center2.x += p[0];
        center2.y += p[1];
        center2.z += p[2];
    });
    center1.divideScalar(inputLandmarksArray.length);
    center2.divideScalar(inputLandmarksArray.length);

    inputLandmarksArray.forEach((value, index) => {
        inputLandmarksArray[index][0] -= center1.x;
        inputLandmarksArray[index][1] -= center1.y;
        inputLandmarksArray[index][2] -= center1.z;
    });

    verticesCopy.forEach((value, index) => {
        verticesCopy[index][0] -= center2.x;
        verticesCopy[index][1] -= center2.y;
        verticesCopy[index][2] -= center2.z;
    });

    // landmarks and vertices are moved to (0,0,0), scale vertices to landmarks.
    // the scale should be 3d vector, but now we assume it can be 1d
    let s1 = 0;
    let s2 = 0;
    inputLandmarksArray.forEach((value) => {
        s1 += value[0] ** 2 + value[1] ** 2 + value[2] ** 2
    });
    s1 = Math.sqrt(s1);

    verticesCopy.forEach((value) => {
        s2 += value[0] ** 2 + value[1] ** 2 + value[2] ** 2
    });
    s2 = Math.sqrt(s2);
    const scale = s1 / s2;
    verticesCopy.forEach((value, index) => {
        verticesCopy[index][0] *= scale;
        verticesCopy[index][1] *= scale;
        verticesCopy[index][2] *= scale;
    });

    // solve rotate
    const landmarksArray = Object.values(inputLandmarksArray);

    const landmarkMat = matrix(landmarksArray);
    const verticesMat = matrix(matrix(verticesCopy).trans());
    const designMat = verticesMat.prod(landmarkMat);

    // TODO:check if rotation need trans()!!!
    const rotation = matrix(computeOptimalRotation(designMat).trans());
    const translate = rotation.prod(matrix([[center2.x * scale], [center2.y * scale], [center2.z * scale]]));

    translate[0][0] = center1.x - translate[0][0];
    translate[1][0] = center1.y - translate[1][0];
    translate[2][0] = center1.z - translate[2][0];

    // console.log("output scale: ", scale);
    // console.log("output r: ", rotation());
    // console.log("output t: ", translate);

    let error = 0;
    const sR = rotation().map((value) => [
        scale * value[0], scale * value[1], scale * value[2]
    ]);
    const inputVerticesMat = matrix(matrix(verticesCopy2).trans());
    const sR_prod_vT = matrix(matrix(sR).prod(inputVerticesMat)).trans();
    sR_prod_vT.forEach((value, index) => {
        sR_prod_vT[index][0] += translate[0][0];
        sR_prod_vT[index][1] += translate[1][0];
        sR_prod_vT[index][2] += translate[2][0];
    });

    sR_prod_vT.forEach((value, index) => {
        error += Math.abs(sR_prod_vT[index][0] - inputLandmarksArrayCopy[index][0]);
        error += Math.abs(sR_prod_vT[index][1] - inputLandmarksArrayCopy[index][1]);
        error += Math.abs(sR_prod_vT[index][2] - inputLandmarksArrayCopy[index][2]);
    });
    // console.log("procrustes error", error);

    return [scale, rotation(), translate]
}

function FrobeniusNorm(mat) {
    return Math.sqrt(mat.reduce((accu0, row) => {
        return accu0 + row.reduce((accu1, ele) => {
            return accu1 + ele * ele;
        }, 0)
    }, 0));
}

function computeOptimalRotation(designMat) {
    if (FrobeniusNorm(designMat) < 1e-9) {
        console.log('Design matrix norm is too small!');
        return null;
    }
    const { u, v, q } = SVD(designMat);
    let postRotation = matrix(u);
    const preRotation = matrix(matrix(v).trans());
    if (postRotation.det() * preRotation.det() < 0) {
        const postArr = postRotation();
        postArr[0][2] *= -1;
        postArr[1][2] *= -1;
        postArr[2][2] *= -1;
        postRotation = matrix(postArr);
    }
    return matrix(postRotation.prod(preRotation));
}

export {
    procrustes
};