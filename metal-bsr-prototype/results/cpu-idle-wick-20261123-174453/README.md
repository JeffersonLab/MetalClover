# One-hour production pause for New Wick contractions

User requested a 3,600-second idle production window. The detached local
controller verifies each process identity before SIGSTOP and before SIGCONT.
It keeps the same native solve, hierarchy, checkpoint and waiting workers.
It automatically resumes after one hour, including if the phone/chat disconnects.
SIGTERM/SIGINT or an explicit matching resume-request.json also releases the
pause. A setup failure resumes every process already stopped in a finally block.

plan.json identifies the exact targets; pause-status.json records the actual
pause time, automatic deadline and eventual duration. The live status is also
/Users/admin/Downloads/MetalClover-results/cfg1020.TIMED_PAUSE.json.
Native solve wall timing includes this idle interval; do not interpret it as
uninterrupted inversion performance. No numerical source, setting or output
is changed, and all standard perambulators remain.
