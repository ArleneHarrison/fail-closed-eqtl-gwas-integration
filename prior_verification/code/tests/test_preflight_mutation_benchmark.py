from run_preflight_mutation_benchmark import base_contract, execute, mutation_registry, run_boundaries


def test_complete_contract_is_ready():
    result, _runtime = execute(base_contract(10))
    assert result["status"] == "READY"
    assert result["status_codes"] == ["OK"]


def test_every_registered_fault_emits_its_expected_code():
    for _name, expected_code, mutation in mutation_registry():
        contract = base_contract(10)
        mutation(contract, 3)
        result, _runtime = execute(contract)
        assert expected_code in result["status_codes"]


def test_prespecified_tolerance_boundaries_match():
    boundary = run_boundaries(10)
    assert boundary["expected_status_matched"].all()
