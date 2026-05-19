"""OpenTelemetry observability for FoldCopilot MCP server.

Graceful no-op when OTel is not installed. Install with:
    pip install foldcopilot[observability]

Follows MCP semantic conventions where applicable:
    https://opentelemetry.io/docs/specs/semconv/gen-ai/mcp/
"""

from __future__ import annotations

import hashlib
import json
import os
import time
from contextlib import contextmanager
from functools import wraps
from typing import Any, Callable, Generator

# --- Graceful OTel import ---

_OTEL_AVAILABLE = False

try:
    from opentelemetry import trace
    from opentelemetry.sdk.trace import TracerProvider
    from opentelemetry.sdk.trace.export import BatchSpanProcessor
    from opentelemetry.exporter.otlp.proto.grpc.trace_exporter import OTLPSpanExporter
    from opentelemetry.sdk.resources import Resource
    from opentelemetry.trace import StatusCode

    _OTEL_AVAILABLE = True
except ImportError:
    trace = None  # type: ignore[assignment]
    StatusCode = None  # type: ignore[assignment,misc]


def is_otel_available() -> bool:
    """Check whether OpenTelemetry SDK is installed."""
    return _OTEL_AVAILABLE


# --- Tracer singleton ---

_tracer: Any = None


def setup_tracing(
    service_name: str = "foldcopilot",
    service_version: str | None = None,
) -> None:
    """Configure the OpenTelemetry TracerProvider with OTLP exporter.

    The OTLP endpoint is read from the standard ``OTEL_EXPORTER_OTLP_ENDPOINT``
    env var (default: ``http://localhost:4317``).

    No-op if OTel is not installed.
    """
    global _tracer

    if not _OTEL_AVAILABLE:
        return

    if service_version is None:
        try:
            import foldcopilot
            service_version = foldcopilot.__version__
        except Exception:
            service_version = "unknown"

    resource = Resource.create(
        {
            "service.name": service_name,
            "service.version": service_version,
        }
    )

    provider = TracerProvider(resource=resource)

    endpoint = os.environ.get("OTEL_EXPORTER_OTLP_ENDPOINT", "http://localhost:4317")
    exporter = OTLPSpanExporter(endpoint=endpoint, insecure=True)
    provider.add_span_processor(BatchSpanProcessor(exporter))

    trace.set_tracer_provider(provider)
    _tracer = trace.get_tracer(service_name, service_version)


def get_tracer() -> Any:
    """Return the configured tracer, or a no-op tracer.

    Safe to call whether or not OTel is installed.
    """
    global _tracer

    if _tracer is not None:
        return _tracer

    if _OTEL_AVAILABLE:
        # Return the global no-op tracer if setup_tracing was never called.
        _tracer = trace.get_tracer("foldcopilot")
        return _tracer

    # OTel not installed — return a lightweight stub.
    return _NoOpTracer()


# --- trace_tool: context manager + decorator ---


def _hash_params(params: dict[str, Any] | None) -> str:
    """Return a short SHA-256 hex digest of the JSON-serialised params."""
    if not params:
        return ""
    try:
        raw = json.dumps(params, sort_keys=True, default=str)
        return hashlib.sha256(raw.encode()).hexdigest()[:16]
    except Exception:
        return ""


@contextmanager
def trace_tool(
    name: str,
    params: dict[str, Any] | None = None,
) -> Generator[Any, None, None]:
    """Context manager that wraps a tool call with an OTel span.

    Attributes set on the span (MCP semantic conventions where applicable):
        - mcp.tool.name
        - mcp.tool.input_params_hash
        - mcp.tool.duration_ms
        - mcp.tool.status  (ok / error)

    Usage::

        with trace_tool("assess_confidence", {"uniprot_id": "P04637"}) as span:
            result = await confidence.assess_confidence(uid)

    Also works as a no-op when OTel is not installed.
    """
    tracer = get_tracer()
    start = time.monotonic()

    if _OTEL_AVAILABLE:
        with tracer.start_as_current_span(f"mcp.tool.{name}") as span:
            span.set_attribute("mcp.tool.name", name)
            params_hash = _hash_params(params)
            if params_hash:
                span.set_attribute("mcp.tool.input_params_hash", params_hash)
            try:
                yield span
                duration_ms = (time.monotonic() - start) * 1000
                span.set_attribute("mcp.tool.duration_ms", round(duration_ms, 2))
                span.set_attribute("mcp.tool.status", "ok")
                span.set_status(StatusCode.OK)
            except Exception as exc:
                duration_ms = (time.monotonic() - start) * 1000
                span.set_attribute("mcp.tool.duration_ms", round(duration_ms, 2))
                span.set_attribute("mcp.tool.status", "error")
                span.set_status(StatusCode.ERROR, str(exc))
                span.record_exception(exc)
                raise
    else:
        # No-op path
        yield _NoOpSpan()


def trace_tool_decorator(name: str | None = None) -> Callable:
    """Decorator form of trace_tool.

    Usage::

        @trace_tool_decorator("assess_confidence")
        async def assess_confidence(uniprot_id: str) -> dict:
            ...
    """

    def decorator(func: Callable) -> Callable:
        tool_name = name or func.__name__

        @wraps(func)
        async def async_wrapper(*args: Any, **kwargs: Any) -> Any:
            with trace_tool(tool_name, kwargs or None):
                return await func(*args, **kwargs)

        @wraps(func)
        def sync_wrapper(*args: Any, **kwargs: Any) -> Any:
            with trace_tool(tool_name, kwargs or None):
                return func(*args, **kwargs)

        import inspect
        if inspect.iscoroutinefunction(func):
            return async_wrapper
        return sync_wrapper

    return decorator


# --- No-op stubs (used when OTel is not installed) ---


class _NoOpSpan:
    """Minimal span stub that silently ignores all calls."""

    def set_attribute(self, key: str, value: Any) -> None:
        pass

    def set_status(self, *args: Any, **kwargs: Any) -> None:
        pass

    def record_exception(self, exc: Exception) -> None:
        pass

    def end(self) -> None:
        pass


class _NoOpTracer:
    """Minimal tracer stub that returns no-op spans."""

    @contextmanager
    def start_as_current_span(self, name: str, **kwargs: Any) -> Generator[_NoOpSpan, None, None]:
        yield _NoOpSpan()

    def start_span(self, name: str, **kwargs: Any) -> _NoOpSpan:
        return _NoOpSpan()
