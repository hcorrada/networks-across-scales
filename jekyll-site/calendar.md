---
layout: default
title: Course Calendar
---

Readings are uploaded to canvas: [https://umd.instructure.com/courses/1030134](https://umd.instructure.com/courses/1030134)

<table class="table">
<colgroup>
<thead>
<tr><th>Date<th>Lecture Materials<th>Readings<th>Work Due
<tbody>
{% for entry in site.data.calendar %}
<tr><td>{{ entry.date }}<td>{% for theme in entry.covers %}[{{ theme.title }}]({{ theme.link }}){% endfor %}<td><td>
{% endfor %}
</table>

Readings are required unless listed *in italics*.

<table class="table">
	<colgroup>
	<thead>
		<tr><th>Legend
	</thead>
	<tbody>
		<tr class="tr-present"><td>Under construction
		<tr class="tr-future"><td>Not updated yet
	</tbody>
</table>
