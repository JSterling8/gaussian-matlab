%% Main function to generate tests
function tests = stage3Tests
    fprintf('Running all tests...\n\n')
    testMandatoryRowSwap();
    testSmallNumbersExampleFromLecture();
    testSuperSmallNumbers();
    testTenByTenZeroValueMatrix();
    testOneByOneMatrix();
    testOneByOneZeroValueMatrix();
    testNonRealNumbersInAFails();
    testNonRealNumbersInbFails()
    testNxMFails();
    testCorrectXValuesReturned100x100();
    testIncorrectSizeOfbFails();
    testPartialRankFails();
    testTwoByTwoMatrix();
    testRandomNByNTest();
    testCorrectXValuesReturned1000x1000();
    fprintf('All tests passed!\n\n')
end

function testOneByOneMatrix()
    fprintf('Beginning 1x1 test.  Will test with 10 random 1x1 matrices\n')
    tic
    
    for i = 1:10
        A = rand(1) .* 100;
        while rank(A) ~= 1
            A = rand(1) .* 100;
        end
        b = rand(1,1) .* 100;
        x = A\b;

        x_calc = stage3(A, b);

        tolerance = 0.00000001;
        for row = 1:1
          if abs(x(row)) - abs (x_calc(row)) > tolerance
              error('Calculated incorrect solution');
          end
        end
    end
    
    toc
    fprintf('1x1 test passed\n\n')
end

function testSmallNumbersExampleFromLecture()
    fprintf('Beginning small numbers test\n')
    tic

    A = [0.001, 0.995; -10.2, 1.00;];
    b = [1; -50];
    x = A\b;
    
    x_calc = stage3(A, b);
    
    tolerance = 0.00000001;
    for row = 1:2
      if abs(x(row)) - abs (x_calc(row)) > tolerance
          error('Calculated incorrect solution');
      end
    end
    
    toc
    fprintf('Small numbers test complete\n\n')
end

function testSuperSmallNumbers()
    fprintf('Beginning super small numbers test\n')
    tic
    A = [0.0000001, 0.00000995; -10.2, 1.00;];
    b = [1; -50];
    x = A\b;
    
    x_calc = stage3(A, b);
    
    tolerance = 0.00000001;
    for row = 1:2
      if abs(x(row)) - abs (x_calc(row)) > tolerance
          error('Calculated incorrect solution');
      end
    end
    
    toc
    fprintf('Super small numbers test complete\n\n')
end

function testOneByOneZeroValueMatrix()
    fprintf('Beginning 1x1 zero value test.\n')
    tic
    
    caught = 0;
    try
    A = zeros(1,1);
    b = zeros(1,1);
    
    x_calc = stage3(A, b);
    catch
        caught = 1;
    end
    
    assert(caught == 1, 'No exception thrown for unsolvable matrix');
    
    toc
    fprintf('1x1 zero value test passed\n\n')
end

function testTenByTenZeroValueMatrix()
    fprintf('Beginning 10x10 zeroes test.\n')
    tic
    
    caught = 0;
    try
    A = zeros(10,10);
    b = zeros(10,1);

    x_calc = stage3(A, b);
    catch
        caught = 1;
    end
    
    assert(caught == 1, 'No exception thrown for unsolvable matrix');
    
    toc
    fprintf('10x10 zeroes test passed\n\n')
end

function testCorrectXValuesReturned100x100()
    fprintf('Beginning 100x100 test.  Will test with 10 random 100x100 matrices\n')
    tic
    
    for i = 1:10
        A = rand(100) .* 100;
        while rank(A) ~= 100
            A = rand(100) .* 100;
        end
        b = rand(100,1) .* 100;
        x = A\b;

        x_calc = stage3(A, b);

        tolerance = 0.00000001;
        for row = 1:100
          if abs(x(row)) - abs (x_calc(row)) > tolerance
              error('Calculated incorrect solution');
          end
        end
    end
    
    toc
    fprintf('100x100 test passed\n\n')
end

function testMandatoryRowSwap()
    fprintf('Beginning mandatory row swap test\n')
    tic
    % The following won't have a solution unless a row swap occurs:
    A = [0, 0, 1; 1, 2, 3; 4, 5, 6];
    % Row one will have to be moved to the bottom.
    b = [4; 3; 2];
    x = A\b;
    
    x_calc = stage3(A, b);
    
    tolerance = 0.00000001;
    for row = 1:3
      if abs(x(row)) - abs (x_calc(row)) > tolerance
          error('Calculated incorrect solution');
      end
    end
    
    toc
    fprintf('Mandatory row swap test complete\n\n')
end

function testPartialRankFails()
    fprintf('Beginning partial rank test\n')
    tic
    
    caught = 0;
    try
        A = [1, 1, 1; 2.5, 2.5, 2.5; 4, 5, 6];
        b = [4; 3; 2];

        stage3(A, b);
    catch
        caught = 1;
    end
    
    assert(caught == 1, 'No error thrown for input with partial rank!');
    
    toc
    fprintf('Partial rank test complete\n\n')
end

function testTwoByTwoMatrix()
    fprintf('Beginning 2x2 test.  Will test with 10 random 2x2 matrices\n')
    tic
    
    for i = 1:10
        A = rand(2) .* 100;
        while rank(A) ~= 2
            A = rand(1) .* 100;
        end
        b = rand(2,1) .* 100;
        x = A\b;

        x_calc = stage3(A, b);

        tolerance = 0.00000001;
        for row = 1:2
          if abs(x(row)) - abs (x_calc(row)) > tolerance
              error('Calculated incorrect solution');
          end
        end
    end
    
    toc
    fprintf('2x2 test passed\n\n')
end

function testNxMFails()
    fprintf('Beginning NxM test (asserting error thrown).\n')
    tic
    
    caught = 0;
    try
        A = zeros(5,4);
        b = zeros(5,1);
        
        stage3(A, b);
    catch
        caught = 1;
    end
    assert(caught == 1, 'No error was thrown for NxM matrix input');
    
    toc
    fprintf('NxM test passed\n\n')
end

function testIncorrectSizeOfbFails()
    fprintf('Beginning b size test (checks both too big and too small) (asserting error thrown).\n')
    tic
    
    caught = 0;
    try
        A = zeros(5,5);
        b = zeros(4,1);
        stage3(A, b);
    catch
        caught = 1;
    end
    assert(caught == 1, 'No error was thrown for invalid b vector input (b too small)');
    
    try
        A = zeros(5,5);
        b = zeros(7,1);
        stage3(A, b);
    catch
        caught = 1;
    end
    assert(caught == 1, 'No error was thrown for invalid b vector input (b too big)');
    
    toc
    fprintf('b size test passed\n\n')
end

function testNonRealNumbersInAFails()
    fprintf('Beginning non-real numbers in A test (asserting error thrown).\n')
    tic
    
    caught = 0;
    try
        A = zeros(5,5);
        A(3,3) = 5i;
        b = zeros(5,1);
        stage3(A, b);
    catch
        caught = 1;
    end
    assert(caught == 1, 'No error was thrown for non-real matrix input');
    
    toc
    fprintf('Non-real numbers in A test passed\n\n')
end

function testNonRealNumbersInbFails()
    fprintf('Beginning non-real numbers in b test (asserting error thrown).\n')
    tic
    
    caught = 0;
    try
        A = zeros(5,5);
        b = zeros(5,1);
        b(1) = 1i;
        stage3(A, b);
    catch
        caught = 1;
    end
    assert(caught == 1, 'No error was thrown for non-real vector input');
    
    toc
    fprintf('Non-real numbers in b test passed\n\n')
end

function testRandomNByNTest()
    fprintf('Beginning random NxN test.  Will test with 100 random NxN matrices\n')
    tic
    
    for i = 1:100
        size = floor(rand(1) * 100);

        A = rand(size) .* 100;
        while rank(A) ~= size
            A = rand(size) .* 100;
        end
        b = rand(size,1) .* 100;
        x = A\b;

        x_calc = stage3(A, b);

        tolerance = 0.00000001;
        for row = 1:size
          if abs(x(row)) - abs (x_calc(row)) > tolerance
              error('Calculated incorrect solution');
          end
        end
    end
    
    toc
    fprintf('Random NxN test passed\n\n')
end

function testCorrectXValuesReturned1000x1000()
    fprintf('Beginning 1000x1000 test.  Will test with 1 random 1000x1000 matrices\n')
    tic
    
    A = rand(1000) .* 100;
    while rank(A) ~= 1000
        A = rand(1000) .* 100;
    end
    b = rand(1000,1) .* 100;
    x = A\b;

    x_calc = stage3(A, b);

    tolerance = 0.00000001;
    for row = 1:1000
      if abs(x(row)) - abs (x_calc(row)) > tolerance
          error('Calculated incorrect solution');
      end
    end
    
    toc
    fprintf('1000x1000 test passed\n\n')
end