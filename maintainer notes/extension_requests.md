## Notes on package extension request:

For us to incorporate a function in the package, we need to be convinced that (1) it is useful to more than a very narrow group of users, and (2) that we thoroughly understand its underpinnings.

For (1), can you describe how broadly such a function might be used?

Point (2) is especially important, because any questions about the use of the function, or any bug reports, come to us, the maintainers of the package. can you give a literature reference for the method you implemented to pool a list of emmGrid objects and obtain their degrees freedom?

Some additional considerations:

- How easily could it be inadvertently misused? I'm thinking for example of a user running emmeans(model, pairwise ~ treatment) and then passing the result to this function.

- What about non-estimability issues? Can this arise? How do we handle this? 

- Generally, what else can go wrong?
