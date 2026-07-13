import { defineCollection } from "astro:content";
import { glob } from "astro/loaders";
import { z } from "astro/zod";

const bioinformatics = defineCollection({
  loader: glob({
    base: "./src/content/bioinformatics",
    pattern: "**/*.{md,mdx}",
  }),

  schema: z.object({
    title: z.string(),
    description: z.string(),
    date: z.coerce.date(),
    updated: z.coerce.date().optional(),
    category: z.string(),
    tags: z.array(z.string()).default([]),
    status: z.string().optional(),
    draft: z.boolean().default(false),
  }),
});

export const collections = {
  bioinformatics,
};