---
name: paper-method-extraction
description: Extract computational biomedical methods, named databases, and software from a research paper. Use when asked to identify reusable analysis tasks from paper text.
---

# Extract Methods from a Biomedical Paper

Use the paper text supplied by the user. Identify computational methods that the paper actually uses and that have specific inputs and outputs. Keep wet-lab procedures, paper-specific hypotheses, and vague labels such as “statistical analysis” out of the task list. Include a database or software package only when the paper names it.

Return one JSON object with exactly these top-level arrays:

- `tasks`: objects with `task_name`, `description`, `inputs`, `outputs`, `code_implementation`, `frequency`, `standard_methods`, and `example`.
- `databases`: objects with `name`, `description`, `url`, `usage`, and `example`.
- `software`: objects with `name`, `description`, `url`, `usage`, and `example`.

For each task, name the actual method, state the data and parameters it needs, state the results it produces, and describe a plausible implementation using established software. Distinguish a suggested implementation from software the paper reports using. In `example`, point to the passage or section showing how this paper used the method. Record `frequency` only when supported by evidence beyond this one paper; otherwise use `null`. For databases and software, give their role in this paper and use `null` for a URL the text does not provide or that has not been verified.

Check every entry against the source text before returning it. Do not turn a proposed analysis into one the authors performed. If the provided text does not show enough detail for an input, output, method, or usage claim, use `null` for that field; omit a task when its method or purpose cannot be identified. Empty categories are `[]`, and the response must remain valid JSON.

Source of the selection criteria and field names: Biomni upstream commit `400c1f3`, `biomni/agent/env_collection.py`, `PaperTaskExtractor.configure()` and `_consolidate_tasks()`. This skill preserves the paper-reading guidance from those prompts; it does not implement or validate the upstream LLM pipeline.
