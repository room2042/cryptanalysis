Hashing Attacks
===============
The :mod:`cryptanalysis.hash` module implements several hashing algorithms that
are vulnerable to a length extension attack. 

Length Extension Attack
-----------------------
A length extension attack on a hash function ``h`` can be used to bypass an
insecure :abbr:`MAC (Message Authentication Code)` scheme. An example of such
an insecure scheme can be described as follows. Generate an authentication code
for a message ``m``, :math:`\sigma_m`, by computing the hash value of a secret
``s`` concatenated with the message ``m``, :math:`\sigma_m =H(s \| m)`. The
code :math:`\sigma_m` for some message ``m`` can be verified using the same
operation. A length extension attack allows us to compute a valid MAC
:math:`\sigma_{m'}` for a message :math:`m' = m \| m_\text{extension}` from a
MAC :math:`\sigma_m` without needing to knowing the secret ``s`` or the message
``m``. However, the length of the originally hashed value, :math:`s \| m` in
this case, needs to be known (or guessed) and :math:`m_\text{extension}` will
be preceded by several bytes that cannot be changed.

To show how a length extension attack works, assume we're given a MAC for a
message ``m = "id=43"`` and we want to obtain a MAC for ``m' = "id=43[arbitrary
bytes]&admin=true"``.

.. testcode:: length-extension-attack

   from cryptanalysis.hash import SHA256

   # obtain the MAC for message m = "id=43", here we assume it's MAC_m
   MAC_m = "6ab7c61b2381e495b3379b025f449b1ed308ea5230377ec217bf8534b356a555"

Using the :meth:`~cryptanalysis.hash.SHA256.extend` method, we obtain the
(candidate) state for the hash object that produces the digest ``MAC_m``, if we
can guess the length of the value that is hashed.

.. testcode:: length-extension-attack

   h = SHA256()
   # we guess that the secret is 16 bytes long, we know that id=43 is 5 bytes
   guessed_length = (16 + 5) * 8
   padding = h.padding(guessed_length)
   m_forge = b"&admin=true"
   for candidate in h.extend(MAC_m, guessed_length):
       candidate.update(m_forge)
       print("m_extension =", padding + m_forge)
       print("forged MAC_m' =", candidate.hexdigest())

   # To verify:
   h = SHA256(b"~SuperSecretKey~id=43")
   print("MAC_m =", h.hexdigest())
   h.update(padding + m_forge)
   print("MAC_m' =", h.hexdigest())

The output of the code shows that our attack was successful.

.. testoutput:: length-extension-attack
   :options: +NORMALIZE_WHITESPACE

   m_extension = bytearray(b'\x80\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\xa8&admin=true')
   forged MAC_m' = 96484ea2b91a14ae7468d673727f9d6c01b0df14d997fd602b6dd4e76ea1a8fc
   MAC_m = 6ab7c61b2381e495b3379b025f449b1ed308ea5230377ec217bf8534b356a555
   MAC_m' = 96484ea2b91a14ae7468d673727f9d6c01b0df14d997fd602b6dd4e76ea1a8fc
