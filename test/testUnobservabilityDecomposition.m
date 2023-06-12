% entry point
function tests = testUnobservabilityDecomposition
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    Ts = 1e-3;
    testCase.TestData.Ts = Ts;
    testCase.TestData.A = eye(5) + Ts * ... 
            [-1.3     1     17      5     16;
               2     -6      3     -1     -8;
               0      1     -8     -7      3;
               2    -13     -2    -15      5;
               6     -6      1      2     -1];
 
    testCase.TestData.B = Ts * [eye(2); zeros(3,2)];

    testCase.TestData.E = Ts * [1   1   1   0   0;
                                0   0   0   1   1]';
 
    testCase.TestData.C = ...
                [1   0   0   0   0;
                 0   1   0   0   0;
                 0   0   1   1   0;
                 0   0   0   0   1];
end

function canConstructTest(testCase)
    td = testCase.TestData;
    sys = LtiProperSystem(td.A, td.B, td.C, td.E);

    testCase.fatalAssertEqual(sys.A, td.A);
    testCase.fatalAssertEqual(sys.B, td.B);
    testCase.fatalAssertEqual(sys.C, td.C);
    testCase.fatalAssertEqual(sys.E, td.E);
    testCase.fatalAssertEqual(size(td.A,1), sys.n);
    testCase.fatalAssertEqual(size(td.C,1), sys.p);
end

function correctInfACinvariantTest(testCase)
    td = testCase.TestData;
    sys = LtiProperSystem(td.A, td.B, td.C, td.E).setInputMap('B');
    W = 1e-3 * [1 0 0 0 0; 0 1 0 0 0]';

    testCase.assertEqual(sys.getInfACinvariantSuspaceContainingInputSpace(), ...
        W, "AbsTol", 1e-6);
end

function correctInfUnobsSpaceTest(testCase)
    td = testCase.TestData;
    sys = LtiProperSystem(td.A, td.B, td.C, td.E).setInputMap('B');
    S = 1e-3 * [1 0 0 0 0; 0 1 0 0 0]';

    testCase.assertEqual(sys.getInfUnobservabilitySubspaceContainingInput(), ...
        S, "AbsTol", 1e-6);
end

function correctTransformTest(testCase)
    td = testCase.TestData;
    sys = LtiProperSystem(td.A, td.B, td.C, td.E).setInputMap('B');

    Texp = [0, 0, 0, 0.0010, 0;
            0, 0, 0, 0,      0.0010;
            1, 0, 0, 0,      0;
            0, 1, 0, 0,      0;
            0, 0, 1, 0,      0];
    Gamexp = [...
             0     0     1     0;
             0     0     0     1;
             1     0     0     0;
             0     1     0     0];
    n1exp = 3;
    p1exp = 2;

    [T,Gamma, n1, p1] = sys.getUnobSubspaceTransform();

    testCase.assertEqual(T, Texp, "AbsTol", 1e-6);
    testCase.assertEqual(Gamma, Gamexp, "AbsTol", 1e-6);
    testCase.assertEqual(n1, n1exp, "AbsTol", 1e-6);
    testCase.assertEqual(p1, p1exp, "AbsTol", 1e-6);
end

function placePolesTest(testCase)
    td = testCase.TestData;
    sys = LtiProperSystem(td.A, td.B, td.C, td.E).setInputMap('B');
    usPoles = -10*[0.1 0.2];
    usModXPoles = -10*[0.5 0.6 0.7];
    L = sys.placePoles(usPoles, usModXPoles);
    exp = sort(reshape([usPoles usModXPoles], [], 1));
    got = sort(reshape(eig(td.A + L*td.C), [], 1));
    testCase.assertEqual(got, exp, "RelTol", 0.1);
end

% function computeGainLQRTest(testCase)
%     td = testCase.TestData;
%     sys = LtiProperSystem(td.A, td.B, td.C, td.E).setInputMap('B');
%     Gexp = [...
%              0         0         0         0;
%              0         0         0         0;
%         0.0000   -0.0010   -0.6695    0.0274;
%        -0.0020    0.0130   -0.3070   -0.0349;
%        -0.0060    0.0060   -0.0013   -0.9892];
%     
%     testCase.assertEqual(sys.designGainLQR(1e2*eye(size(td.A, 1)), ...
%                     eye(size(td.C, 1))), Gexp, "AbsTol", 5e-5);
% end

% function testBlocks(testCase)
%     td = testCase.TestData;
%     sys = LtiProperSystem(td.A, td.B, td.C, td.E).setInputMap('B');
%     sys.designGainLQR(1e2*eye(size(td.A, 1)), eye(size(td.C, 1)))
% 
%     A11 = [ 0.3159   -0.6831    0.0306;
%            -0.3121    0.6749   -0.0302;
%            -0.0003    0.0007   -0.0000];
%     A21 = [ 17     5    16;
%              3    -1    -8];
%     A22 = [ 0.9987    0.0010
%             0.0020    0.9940];
%     B2 = eye(2);
%     C11 = [[1;0] eye(2)];
%     C21 = zeros(2,3);
%     C22 = 1e-3 * eye(2);
%     E1 = [eye(2); [0 1]] * 1e-3;
%     E2 = [1 0; 1 0];
% 
%     syst = sys.getBlockUnobSubspaceDecomposition();
% 
%     testCase.assertEqual(syst.A.getBlock(1,1), A11, "AbsTol", 5e-5);
%     testCase.assertEqual(syst.A.getBlock(2,1), A21, "AbsTol", 5e-5);
%     testCase.assertEqual(syst.A.getBlock(2,2), A22, "AbsTol", 5e-5);
%     testCase.assertEqual(syst.B.getBlock(2,1),  B2, "AbsTol", 5e-5);
%     testCase.assertEqual(syst.C.getBlock(1,1), C11, "AbsTol", 5e-5);
%     testCase.assertEqual(syst.C.getBlock(2,1), C21, "AbsTol", 5e-5);
%     testCase.assertEqual(syst.C.getBlock(2,2), C22, "AbsTol", 5e-5);
%     testCase.assertEqual(syst.E.getBlock(1,1),  E1, "AbsTol", 5e-5);
%     testCase.assertEqual(syst.E.getBlock(2,1),  E2, "AbsTol", 5e-5);
% end
