"""Test that seq_ids are unique when direction='both'.

Verifies that step0_discover_dmrs adds direction suffixes (_U, _M) to seq_ids
when running both hypo and hyper directions, preventing duplicate seq_ids
that cause wrong DMR metadata lookups in Excel output.
"""
import json
import os
import tempfile
from unittest.mock import patch, MagicMock
from collections import Counter


def make_mock_block(cell_type_id, rank, chrom, start, end, direction):
    """Create a mock DMRBlock-like object with the given fields."""
    block = MagicMock()
    block.cell_type_id = cell_type_id
    block.seq_id = f"{cell_type_id}_{rank:04d}"
    block.rank = rank
    block.chrom = chrom
    block.start = start
    block.end = end
    block.num_cpgs = 5
    block.block_len = end - start
    block.gene = ""
    block.annotation = ""
    block.dms_score = 0.7
    block.target_mean_meth = 0.9
    block.background_mean_meth = 0.1
    block.delta_means = 0.8
    block.delta_quants = 0.7
    block.delta_maxmin = 0.85
    block.ttest_pval = 0.001
    block.mwtest_pval = 0.002
    block.mvalue_ttest = 0.001
    block.direction = direction
    block.target_celltype = cell_type_id
    block.bg_subgroup_meth = {}
    block.cpg_sites = []
    return block


class TestSeqIdUniqueness:
    """When direction='both', hypo and hyper runs produce blocks with the same
    ranks, so seq_ids would collide without the direction suffix fix."""

    def test_no_duplicate_seq_ids_with_both_directions(self):
        """step0_discover_dmrs with direction='both' should produce unique seq_ids."""
        from methyl_panel.pipeline import step0_discover_dmrs

        # Mock blocks: hypo and hyper both have ranks 1-3
        hypo_blocks = [
            make_mock_block("Blood-T-CD3", 1, "chr1", 1000, 1200, "U"),
            make_mock_block("Blood-T-CD3", 2, "chr2", 2000, 2200, "U"),
            make_mock_block("Blood-T-CD3", 3, "chr3", 3000, 3200, "U"),
        ]
        hyper_blocks = [
            make_mock_block("Blood-T-CD3", 1, "chr4", 4000, 4200, "M"),
            make_mock_block("Blood-T-CD3", 2, "chr5", 5000, 5200, "M"),
            make_mock_block("Blood-T-CD3", 3, "chr6", 6000, 6200, "M"),
        ]

        args = MagicMock()
        args.cell_type = "Blood-T-CD3"
        args.direction = "both"
        args.beta_dir = "/fake/beta"
        args.blocks_file = None
        args.groups_csv = None
        args.output_dir = tempfile.mkdtemp()
        args.threads = 1
        args.top_markers = 100
        args.delta_means = 0.3
        args.unmeth_mean_thresh = 0.2
        args.meth_mean_thresh = 0.5
        args.min_bg_subgroup_meth = 0.0
        args.wgbstools_path = None
        args.max_bg_samples = None
        args.skip_find_markers = True
        args.use_full_atlas = True

        # Mock discover_dmrs to return our fake blocks
        call_count = [0]
        def mock_discover_dmrs(**kwargs):
            call_count[0] += 1
            return hypo_blocks if call_count[0] == 1 else hyper_blocks

        with patch("methyl_panel.phase0_dmr_discovery.discover_dmrs",
                   side_effect=mock_discover_dmrs):
            with patch("methyl_panel.phase0_dmr_discovery.VALID_CELL_TYPES",
                       ["Blood-T-CD3"]):
                result = step0_discover_dmrs(args)

        # Check that all seq_ids are unique
        seq_ids = [b["seq_id"] for b in result]
        id_counts = Counter(seq_ids)
        dups = {k: v for k, v in id_counts.items() if v > 1}

        assert len(dups) == 0, f"Duplicate seq_ids found: {dups}"
        assert len(seq_ids) == 6, f"Expected 6 blocks, got {len(seq_ids)}"

        # Verify suffixes are correct
        hypo_ids = [b["seq_id"] for b in result if b.get("direction") == "U"]
        hyper_ids = [b["seq_id"] for b in result if b.get("direction") == "M"]

        for sid in hypo_ids:
            assert sid.endswith("_U"), f"Hypo seq_id should end with _U: {sid}"
        for sid in hyper_ids:
            assert sid.endswith("_M"), f"Hyper seq_id should end with _M: {sid}"

    def test_no_suffix_for_single_direction(self):
        """When direction is not 'both', seq_ids should NOT have suffixes."""
        from methyl_panel.pipeline import step0_discover_dmrs

        hypo_blocks = [
            make_mock_block("Blood-T-CD3", 1, "chr1", 1000, 1200, "U"),
            make_mock_block("Blood-T-CD3", 2, "chr2", 2000, 2200, "U"),
        ]

        args = MagicMock()
        args.cell_type = "Blood-T-CD3"
        args.direction = "hypo"
        args.beta_dir = "/fake/beta"
        args.blocks_file = None
        args.groups_csv = None
        args.output_dir = tempfile.mkdtemp()
        args.threads = 1
        args.top_markers = 100
        args.delta_means = 0.3
        args.unmeth_mean_thresh = 0.2
        args.meth_mean_thresh = 0.5
        args.min_bg_subgroup_meth = 0.0
        args.wgbstools_path = None
        args.max_bg_samples = None
        args.skip_find_markers = True
        args.use_full_atlas = True

        with patch("methyl_panel.phase0_dmr_discovery.discover_dmrs",
                   return_value=hypo_blocks):
            with patch("methyl_panel.phase0_dmr_discovery.VALID_CELL_TYPES",
                       ["Blood-T-CD3"]):
                result = step0_discover_dmrs(args)

        seq_ids = [b["seq_id"] for b in result]
        for sid in seq_ids:
            assert not sid.endswith("_U"), (
                f"Single-direction seq_id should NOT have _U suffix: {sid}"
            )
            assert not sid.endswith("_M"), (
                f"Single-direction seq_id should NOT have _M suffix: {sid}"
            )


class TestSeqIdValidationGuard:
    """The _validate_seq_id_uniqueness guard prevents stale dmr_blocks.json
    files (generated before the v2.2.21 fix) from being used in downstream
    steps."""

    def test_validation_raises_on_duplicate_seq_ids(self):
        """_validate_seq_id_uniqueness should raise RuntimeError on duplicates."""
        from methyl_panel.pipeline import _validate_seq_id_uniqueness

        blocks = [
            {"seq_id": "CellType_0001", "chrom": "chr1"},
            {"seq_id": "CellType_0001", "chrom": "chr2"},  # duplicate
            {"seq_id": "CellType_0002", "chrom": "chr3"},
        ]

        try:
            _validate_seq_id_uniqueness(blocks, "/fake/path/dmr_blocks.json")
            assert False, "Should have raised RuntimeError"
        except RuntimeError as e:
            assert "CellType_0001" in str(e)
            assert "v2.2.21" in str(e)
            assert "--discover-dmrs" in str(e)

    def test_validation_passes_on_unique_seq_ids(self):
        """_validate_seq_id_uniqueness should pass when all seq_ids are unique."""
        from methyl_panel.pipeline import _validate_seq_id_uniqueness

        blocks = [
            {"seq_id": "CellType_0001_U", "chrom": "chr1"},
            {"seq_id": "CellType_0001_M", "chrom": "chr2"},
            {"seq_id": "CellType_0002_U", "chrom": "chr3"},
        ]

        # Should not raise
        _validate_seq_id_uniqueness(blocks)

    def test_step2_rejects_stale_dmr_blocks(self):
        """step2_bisulfite_convert should reject a stale dmr_blocks.json
        with duplicate seq_ids when loading from disk."""
        from methyl_panel.pipeline import step2_bisulfite_convert

        # Create a stale dmr_blocks.json with duplicate seq_ids
        tmp_dir = tempfile.mkdtemp()
        stale_blocks = [
            {"seq_id": "CT_0001", "chrom": "chr1", "start": 1000, "end": 1200},
            {"seq_id": "CT_0001", "chrom": "chr2", "start": 2000, "end": 2200},
        ]
        blocks_path = os.path.join(tmp_dir, "dmr_blocks.json")
        with open(blocks_path, 'w') as f:
            json.dump(stale_blocks, f)

        args = MagicMock()
        args.output_dir = tmp_dir
        args.max_blocks = None
        args.flank = 200

        try:
            step2_bisulfite_convert(args, blocks_data=None)
            assert False, "Should have raised RuntimeError"
        except RuntimeError as e:
            assert "Duplicate seq_ids" in str(e)
        finally:
            import shutil
            shutil.rmtree(tmp_dir)
