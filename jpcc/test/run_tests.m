function run_tests()
try
    test_conservation_3step_alt()
catch ME
    rethrow(ME);
end
end