"""Tests for TMalphaFold client — membrane orientation data."""

from unittest.mock import AsyncMock, patch

import pytest

from foldcopilot.clients.tmalphaFold_client import (
    get_membrane_context,
    get_membrane_topology,
)


class TestGetMembraneTopology:
    @pytest.mark.asyncio
    async def test_not_found(self):
        """Non-TM protein returns found=False."""
        import httpx
        from unittest.mock import MagicMock

        mock_request = MagicMock()
        mock_response = AsyncMock()
        mock_response.status_code = 404

        def raise_for_status():
            raise httpx.HTTPStatusError("404", request=mock_request, response=mock_response)

        mock_response.raise_for_status = raise_for_status

        with patch("foldcopilot.clients.tmalphaFold_client.httpx.AsyncClient") as mock_client:
            mock_instance = AsyncMock()
            mock_instance.get.return_value = mock_response
            mock_instance.__aenter__ = AsyncMock(return_value=mock_instance)
            mock_instance.__aexit__ = AsyncMock(return_value=False)
            mock_client.return_value = mock_instance

            # Clear cache
            from foldcopilot.clients.tmalphaFold_client import _cache_path
            cache = _cache_path("P04637", "topology")
            if cache.exists():
                cache.unlink()

            result = await get_membrane_topology("P04637")
            assert result["found"] is False

    @pytest.mark.asyncio
    async def test_found(self):
        """TM protein returns segments."""
        from unittest.mock import MagicMock

        mock_response = MagicMock()  # httpx Response methods are all sync
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "type": "alpha-helical",
            "tm_segments": [
                {"start": 30, "end": 55, "type": "transmembrane"},
                {"start": 70, "end": 95, "type": "transmembrane"},
            ],
            "tilt_angle": 15.3,
        }

        with patch("foldcopilot.clients.tmalphaFold_client.httpx.AsyncClient") as mock_client:
            mock_instance = AsyncMock()
            mock_instance.get.return_value = mock_response
            mock_instance.__aenter__ = AsyncMock(return_value=mock_instance)
            mock_instance.__aexit__ = AsyncMock(return_value=False)
            mock_client.return_value = mock_instance

            # Clear any cache
            from foldcopilot.clients.tmalphaFold_client import _cache_path
            cache = _cache_path("P08172", "topology")
            if cache.exists():
                cache.unlink()

            result = await get_membrane_topology("P08172")
            assert result["found"] is True
            assert result["n_tm_segments"] == 2
            assert result["topology_type"] == "alpha-helical"
            assert result["membrane_insertion_angle"] == 15.3


class TestGetMembraneContext:
    @pytest.mark.asyncio
    async def test_returns_combined(self):
        """Combined context has both TMalphaFold and OPM."""
        with patch("foldcopilot.clients.tmalphaFold_client.get_membrane_topology") as mock_tm, \
             patch("foldcopilot.clients.tmalphaFold_client.get_opm_orientation") as mock_opm:
            mock_tm.return_value = {
                "found": True,
                "n_tm_segments": 7,
                "topology_type": "alpha-helical",
            }
            mock_opm.return_value = {"found": False}

            result = await get_membrane_context("P08172")
            assert result["is_transmembrane"] is True
            assert "tmalphaFold" in result
            assert "opm" in result
            assert "7 transmembrane" in result["interpretation"]

    @pytest.mark.asyncio
    async def test_not_tm(self):
        """Soluble protein returns is_transmembrane=False."""
        with patch("foldcopilot.clients.tmalphaFold_client.get_membrane_topology") as mock_tm, \
             patch("foldcopilot.clients.tmalphaFold_client.get_opm_orientation") as mock_opm:
            mock_tm.return_value = {"found": False, "n_tm_segments": 0}
            mock_opm.return_value = {"found": False}

            result = await get_membrane_context("P04637")
            assert result["is_transmembrane"] is False
            assert "soluble" in result["interpretation"].lower()
