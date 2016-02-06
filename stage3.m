function [ x ] = stage3( A, b )
%stage1 Solves for x in equation Ax = b
%   Returns a solution vector given a coefficient matrix and
%   right-hand-side values

    % Check that A is square
    dimensions = size(A);
    if dimensions(1,1) ~= dimensions(1,2)
        error('Input coefficient matrix is not square')
    end

    % Check that b is the right size
    if length(b) ~= dimensions(1,1)
       error('Input values vector is not the correct size') 
    end

    % Check that A contains only real numbers
    for index = 1:numel(A)
        if A(index)*A(index) < 0
            error('Input coefficient matrix has imaginary number(s)')
        end
    end

    % Check that b contains only real numbers
    for index = 1:length(b)
       if b(index)*b(index) < 0
           error('Input values vector has imaginary number(s)');
       end
    end

    row_count = dimensions(1,1);
    column_count = dimensions(1,2);
    AUG = [A,b];

    % For each column, 1->n
    for column_inspecting = 1:column_count
        % Initialize mutator row number to the top possible column...
        mutator_row_number = column_inspecting;

        % Pick the largest pivot/mutator at or below the existing 
        % mutator_row_number
        for row_inspecting = column_inspecting:row_count
            if AUG(row_inspecting, column_inspecting) > AUG(mutator_row_number, column_inspecting)
               mutator_row_number = row_inspecting;
            end
        end

        % Swap the mutator row so the pivot is on the diagonal
        if mutator_row_number == column_inspecting
            % Do nothing... It's already on the diagonal
        else 
           temp_A = AUG(mutator_row_number, :);

           % Move diagonal row to old pivotal row
           AUG(mutator_row_number, :) = AUG(column_inspecting, :);

           % Move pivotal row to diagonal row
           AUG(column_inspecting, :) = temp_A;
        end

        % Set the mutator_row_number to be the diagonal again
        mutator_row_number = column_inspecting;

        % For each row from 2->n
        for row_inspecting = 2:row_count     
            % If the column we're looking at is under the diagonal
            if column_inspecting < row_inspecting
                % Make that cell 0 using the mutator row.  As long as we use
                % a row above the current row, we'll never unset a 0
                % in a previous column
                cell_value_in_current_row = AUG(row_inspecting, column_inspecting);

                % If it's already 0, we don't have to do anything
                if cell_value_in_current_row == 0
                    continue
                end

                cell_value_in_mutator_row = AUG(mutator_row_number, column_inspecting);

                % We can't divide by 0... This state *should* never happen
                if cell_value_in_mutator_row == 0
                    error('Not enough info to convert to upper echelon form.')
                end

                % Find out what we have to divide our mutator row by in order
                % to create a 0 in the row/column inspecting
                multiplication_factor = cell_value_in_current_row / cell_value_in_mutator_row;

                % Create the transformed mutator rows, using our multiplication
                % factor
                mutator_row_A = AUG(mutator_row_number, :) .* multiplication_factor;

                % Subtract our mutator row from the row we're inspecting
                AUG(row_inspecting, :) = AUG(row_inspecting, :) - mutator_row_A;
            end
        end
    end

    % Determine if matrix is full rank by calculating it's determinant (the
    % product of its diagonal)
    diagonal_product = 1;
    for index = 1:row_count
        diagonal_product = diagonal_product * AUG(index, index);
    end
    
    if diagonal_product == 0
        error('Cannot solve.  Matrix does not have full rank')
    else
        fprintf('Matrix has full rank. Solving via back substitution.\n')
    end

    % Solve from last row back to first...
    x = zeros(row_count,1);
    row_solving = row_count;

    while row_solving >= 1
        % If we're solving for X1 in (#1*X1 + #2*X2 + #3*X3 = 4),
        % calculate #2*X2 and #3*X3, subtract that value from both sides of
        % the equation, then divide both sides by #1 to get the value of X1
        
        value_of_row_after_unknown = 0;
        column_index = row_solving + 1;
        
        while column_index <= column_count
            value_of_row_after_unknown = value_of_row_after_unknown + (AUG(row_solving, column_index) * x(column_index));

            column_index = column_index + 1;
        end

        x(row_solving) = (AUG(row_solving, column_count + 1) - value_of_row_after_unknown) / AUG(row_solving, row_solving);

        row_solving = row_solving - 1;
    end
end