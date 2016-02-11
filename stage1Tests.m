
%% Main function to generate tests
function tests = stage1Tests
    fprintf('Running all tests...\n\n')
    
    testOneByOneMatrix();
    testOneByOneZeroValueMatrix();
    testTwoByTwoMatrix();
    test100x100();
    test10x10();
    testTenByTenZeroValueMatrix();
    testNxMFails();
    testNonRealNumbersInAFails();
    testSuperSmallNumbers();
    testMandatoryRowSwap();
    testSparseMatrix();
    testEmptyFails;
    testRandomNByNTest();
    test1000x1000();
   
    fprintf('All tests passed!\n\n')
end

function testOneByOneMatrix()
    fprintf('Beginning 1x1 test.  Will test with 10 random 1x1 matrices\n')
    tic
    
    for i = 1:10
        value = rand(1) .* 100;
        A = value;
        upper = stage1(A);
        
        assert(upper(1) == value);
    end
    
    toc
    fprintf('1x1 test passed\n\n')
end

function testOneByOneZeroValueMatrix()
    fprintf('Beginning 1x1 zero value test.\n')
    tic
    
    A = zeros(1,1);
    
    upper = stage1(A);
    
    assert(upper(1) == 0);
    dimensions = size(upper);
    
    assert(dimensions(1,1) == 1);
    assert(dimensions(1,2) == 1);

    toc
    fprintf('1x1 zero value test passed\n\n')
end

function testTenByTenZeroValueMatrix()
    fprintf('Beginning 10x10 zeroes test.\n')
    tic
    
    A = zeros(10,10);

    upper = stage1(A);

    for index = 1:numel(upper)
        assert(upper(index) == 0, 'Invalid upper output');
    end
    
    toc
    fprintf('10x10 zeroes test passed\n\n')
end

function test10x10()
    fprintf('Beginning 10x00 test.  Will test with 100 random 10x10 matrices\n')
    tic
    
    for i = 1:100
        A = rand(10) .* 100;
        while rank(A) ~= 10 || cond(A) > 10^4
            A = rand(10) .* 100;
        end

        upper = stage1(A);

        % Assert all elements below the diagonal are 0
        for row_index = 2:10
           for column_index = 1:row_index-1
              assert(abs(upper(row_index, column_index)) < 0.000000001); 
           end
        end
    end
    
    toc
    fprintf('10x10 test passed\n\n')
end

function testSparseMatrix()
    fprintf('Beginning sparse matrix test.\n')
    tic

    for index = 1:10
        A = zeros(10,10);
        A(index,1) = index;
        upper = stage1(A);
        
        assert(upper(1,1) == index);
    end
    
    % This also checks for pivoting, ensuring the largest magnitude is on
    % the diagonal
    for index = 1:9
        A = zeros(10,10);
        A(index + 1, index) = index;
        upper = stage1(A);
        
        assert(upper(index, index) == index);
    end
    
    toc
    fprintf('Sparse matrix test passed\n\n');
end

function test100x100()
    fprintf('Beginning 100x100 test.  Will test with 10 random 100x100 matrices\n')
    tic
    
    for i = 1:10
        A = rand(100) .* 100;
        while rank(A) ~= 100 || cond(A) > 10^4
            A = rand(100) .* 100;
        end

        upper = stage1(A);

        % Assert all elements below the diagonal are 0
        for row_index = 2:100
           for column_index = 1:row_index-1
              assert(abs(upper(row_index, column_index)) < 0.000000001); 
           end
        end
    end
    
    toc
    fprintf('100x100 test passed\n\n')
end

function testTwoByTwoMatrix()
    fprintf('Beginning 2x2 test.  Will test with 10 random 2x2 matrices\n')
    tic
    
    for i = 1:10
        A = rand(2) .* 100;
        while rank(A) ~= 2 || cond(A) > 10^4
            A = rand(1) .* 100;
        end

        upper = stage1(A);

        assert(abs(upper(2,1)) < 0.00000001);
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
        
        stage1(A);
    catch
        caught = 1;
    end
    
    assert(caught == 1, 'No error was thrown for NxM matrix input');
    
    toc
    fprintf('NxM test passed\n\n')
end

function testEmptyFails()
    fprintf('Beginning empty input test (asserting error thrown).\n')
    tic
    
    caught = 0;
    try
        A = zeros(0,0);
        
        stage1(A);
    catch
        caught = 1;
    end
    
    assert(caught == 1, 'No error was thrown for empty input');
    
    toc
    fprintf('NxM test passed\n\n')
end

function testSuperSmallNumbers()
    fprintf('Beginning super small numbers test 2\n')
    tic

    A = [0.0000001, 0.00000995; -10.2, 1.00;];
    upper = stage1(A);
    assert(upper(2,1) == 0);
    
    toc
    fprintf('Super small numbers test 2 complete\n\n')
end

function testNonRealNumbersInAFails()
    fprintf('Beginning non-real numbers in A test (asserting error thrown).\n')
    tic
    
    caught = 0;
    try
        A = zeros(5,5);
        A(3,3) = 5i;
        stage1(A);
    catch
        caught = 1;
    end
    assert(caught == 1, 'No error was thrown for non-real matrix input');
    
    toc
    fprintf('Non-real numbers in A test passed\n\n')
end

function testRandomNByNTest()
    fprintf('Beginning random NxN test.  Will test with 100 random NxN matrices\n')
    tic
    
    for i = 1:100
        size = floor(rand(1) * 100) + 1;

        A = rand(size) .* 100;
        while rank(A) ~= size || cond(A) > 10^4
            A = rand(size) .* 100;
        end

        upper = stage1(A);

        
         % Assert all elements below the diagonal are 0
        for row_index = 2:size
           for column_index = 1:row_index-1
              assert(upper(row_index, column_index) < 0.00000001); 
           end
        end
    end
    
    toc
    fprintf('Random NxN test passed\n\n')
end

function test1000x1000()
    fprintf('Beginning 1000x1000 test.  Will test with 1 random 1000x1000 matrices\n')
    tic
    
    A = rand(1000) .* 100;
    while rank(A) ~= 1000 || cond(A) > 10^5
        A = rand(1000) .* 100;
    end
    b = rand(1000,1) .* 100;
    x = A\b;

    upper = stage1(A);

    % Assert all elements below the diagonal are 0
    for row_index = 2:999
       for column_index = 1:row_index-1
          assert(upper(row_index, column_index) < 0.00000001); 
       end
    end
    
    toc
    fprintf('1000x1000 test passed\n\n')
end

function testMandatoryRowSwap()
    fprintf('Beginning mandatory row swap test 2\n')
    tic
    % The following won't have a solution unless a row swap occurs:
    A = [1, 2, 3; 4, 1, 5; 3, 6, 9];
    
    upper = stage1(A);
    
     % Assert all elements below the diagonal are 0
    for row_index = 2:3
       for column_index = 1:row_index-1
          assert(upper(row_index, column_index) == 0); 
       end
    end
    
    toc
    fprintf('Mandatory row swap test 2 complete\n\n')
end