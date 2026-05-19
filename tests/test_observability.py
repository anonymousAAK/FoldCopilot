"""Tests for foldcopilot.observability — graceful no-op behaviour."""

from __future__ import annotations

import importlib
import sys
from unittest import mock

import pytest


# ---------------------------------------------------------------------------
# 1. Module loads without OTel installed (graceful no-op)
# ---------------------------------------------------------------------------


def test_module_loads_without_otel():
    """observability module must import cleanly even when OTel is absent."""
    # Temporarily hide OTel packages so the import falls back to no-op.
    otel_modules = [key for key in sys.modules if key.startswith("opentelemetry")]
    saved = {key: sys.modules.pop(key) for key in otel_modules}

    with mock.patch.dict(sys.modules, {
        "opentelemetry": None,
        "opentelemetry.trace": None,
        "opentelemetry.sdk.trace": None,
        "opentelemetry.sdk.trace.export": None,
        "opentelemetry.exporter.otlp.proto.grpc.trace_exporter": None,
        "opentelemetry.sdk.resources": None,
    }):
        # Force re-import
        if "foldcopilot.observability" in sys.modules:
            del sys.modules["foldcopilot.observability"]
        import foldcopilot.observability as obs
        obs = importlib.reload(obs)

        assert obs.is_otel_available() is False

    # Restore original modules
    sys.modules.update(saved)
    if "foldcopilot.observability" in sys.modules:
        del sys.modules["foldcopilot.observability"]


# ---------------------------------------------------------------------------
# 2. trace_tool works as a no-op (no crash, no side effects)
# ---------------------------------------------------------------------------


def test_trace_tool_noop():
    """trace_tool context manager must work as a no-op without OTel."""
    from foldcopilot.observability import trace_tool, _NoOpSpan

    # Even if OTel IS installed, the context manager should still work.
    with trace_tool("test_tool", {"key": "value"}) as span:
        # span is either a real OTel span or a _NoOpSpan — both are fine
        assert span is not None


def test_trace_tool_captures_exception():
    """trace_tool must re-raise exceptions after recording them."""
    from foldcopilot.observability import trace_tool

    with pytest.raises(ValueError, match="boom"):
        with trace_tool("failing_tool"):
            raise ValueError("boom")


def test_trace_tool_decorator_sync():
    """trace_tool_decorator must work on sync functions."""
    from foldcopilot.observability import trace_tool_decorator

    @trace_tool_decorator("my_sync_tool")
    def my_tool(x: int) -> int:
        return x * 2

    assert my_tool(5) == 10


@pytest.mark.asyncio
async def test_trace_tool_decorator_async():
    """trace_tool_decorator must work on async functions."""
    from foldcopilot.observability import trace_tool_decorator

    @trace_tool_decorator("my_async_tool")
    async def my_tool(x: int) -> int:
        return x * 2

    assert await my_tool(5) == 10


# ---------------------------------------------------------------------------
# 3. setup_tracing doesn't crash without OTel
# ---------------------------------------------------------------------------


def test_setup_tracing_noop_without_otel():
    """setup_tracing must be a silent no-op when OTel is not installed."""
    otel_modules = [key for key in sys.modules if key.startswith("opentelemetry")]
    saved = {key: sys.modules.pop(key) for key in otel_modules}

    with mock.patch.dict(sys.modules, {
        "opentelemetry": None,
        "opentelemetry.trace": None,
        "opentelemetry.sdk.trace": None,
        "opentelemetry.sdk.trace.export": None,
        "opentelemetry.exporter.otlp.proto.grpc.trace_exporter": None,
        "opentelemetry.sdk.resources": None,
    }):
        if "foldcopilot.observability" in sys.modules:
            del sys.modules["foldcopilot.observability"]
        import foldcopilot.observability as obs
        obs = importlib.reload(obs)

        # Must not raise
        obs.setup_tracing()
        assert obs.is_otel_available() is False

        # get_tracer should return a no-op tracer
        tracer = obs.get_tracer()
        assert tracer is not None

    sys.modules.update(saved)
    if "foldcopilot.observability" in sys.modules:
        del sys.modules["foldcopilot.observability"]


def test_setup_tracing_does_not_crash_with_otel():
    """setup_tracing must not crash when OTel IS installed (if present)."""
    from foldcopilot.observability import is_otel_available, setup_tracing

    # This test only validates no crash — doesn't require OTel installed.
    # If OTel is installed, it configures a real provider; if not, it no-ops.
    setup_tracing(service_name="foldcopilot-test")


def test_get_tracer_returns_usable_tracer():
    """get_tracer must return something with start_as_current_span."""
    from foldcopilot.observability import get_tracer

    tracer = get_tracer()
    assert tracer is not None
    assert hasattr(tracer, "start_as_current_span")
