%% Main function to generate tests
function tests = stage3Tests
    testMandatoryRowSwap();
    testTenByTenZeroValueMatrix();
    testOneByOneMatrix();
    testOneByOneZeroValueMatrix();
    testNonRealNumbersFail();
    testNxMFails();
    testCorrectXValuesReturned100x100();
    testCorrectXValuesReturned1000x1000();
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

function testNxMFails()
    fprintf('Beginning NxM test (asserting error thrown).\n')
    tic
    
    caught = 0;
    try
        A = zeros(5,4);
        stage3(A);
    catch
        caught = 1;
    end
    assert(caught == 1, 'No error was thrown for NxM matrix input');
    
    toc
    fprintf('NxM test passed\n\n')
end

function testNonRealNumbersFail()
    fprintf('Beginning non-real number test (asserting error thrown).\n')
    tic
    
    caught = 0;
    try
        A = zeros(5,5);
        A(3,3) = 5i;
        stage3(A);
    catch
        caught = 1;
    end
    assert(caught == 1, 'No error was thrown for non-real matrix input');
    
    toc
    fprintf('Non-real number test test passed\n\n')
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